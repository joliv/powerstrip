#include <stdio.h>   // for printf
#include <cmath>     // for abs, log2
#include <cstdint>   // for fixed ints like uint16_t
#include <fstream>   // for ifstream
#include <utility>   // I don't remember
#include <vector>    // for vector
#include <algorithm> // for max
#include <string>    // for std::string
#include <chrono>

#include "emmintrin.h" // for __m128, intel intrinsic
#include "vendor/simdcomp/include/simdcomp.h" // for bit packing

#include "compress.hpp"

#ifdef DEBUG
#define dbg(...) do {\
    printf("  ");\
    printf(__VA_ARGS__);\
    printf("\n");\
} while(0)
#else
#define dbg(...) do {} while(0);
#endif

#define writebig(buf, offset, size, x) do {\
    memcpy(outbuf + offset, x, size);\
    offset += size;\
} while(0)

// One pass to find the noise floor, basically just the mode
uint16_t get_zeroval(const uint16_t* inbuf, const size_t numints) {
    uint16_t common = 0;
    uint16_t common_count = 0;
    size_t* histogram = (size_t*)calloc(UINT16_MAX+1, sizeof(size_t));
    for (size_t i = 0; i < numints; i++) {
        histogram[inbuf[i]]++;
        if (histogram[inbuf[i]] > common_count) {
            common = inbuf[i];
            common_count = histogram[inbuf[i]];
        }
    }
    free(histogram);
    if (10 * common_count < numints) {
        return common;
    }
    return UINT16_MAX;
}

// Simple delta encoding. See how lemire/simdcomp does this fast?
// (Nevermind: the answer is SIMD, of course)
void delta_encode(int32_t* xs, uint32_t xs_len) {
    int32_t prev = 0;
    for (int i = 0; i < xs_len; i++) {
        int32_t tmpprev = xs[i];
        xs[i] -= prev;
        prev = tmpprev;
    }
}

struct header {
    uint32_t indexsize;
    uint32_t numints;
    uint16_t zeroval;
    uint32_t outliersize;
    uint32_t siglen;
    uint8_t  bitsneeded;
};

#define writesmall(buf, offset, type, x) do {\
    type* interp = reinterpret_cast<type*>(buf+offset);\
    interp[0] = x;\
    offset += sizeof(type);\
} while(0)

size_t writeheader(char* outbuf, size_t offset, struct header h) {
    writesmall(outbuf, offset, uint32_t, h.indexsize);
    writesmall(outbuf, offset, uint32_t, h.numints);
    writesmall(outbuf, offset, uint16_t, h.zeroval);
    writesmall(outbuf, offset, uint32_t, h.outliersize);
    writesmall(outbuf, offset, uint32_t, h.siglen);
    writesmall(outbuf, offset, uint8_t,  h.bitsneeded);
    return offset;
}

int64_t deltacode(uint16_t* inbuf, size_t insize, char* outbuf) {
    const size_t numints = insize / 2; // 2 bytes per uint16
    dbg("Encoding %zu ints", numints);

    // const uint16_t zeroval = 0, zerothresh = 0;
    uint16_t zeroval = get_zeroval(inbuf, numints);
    uint16_t zerothresh = zeroval + 5; // Sure, let's say +5W for now

    // No floor exists.
    if (zeroval == UINT16_MAX) {
        zeroval = 0;
        zerothresh = 0;
        dbg("No floor found. Falling back to lossless");
    } else {
        dbg("Floor found at %d watts", zeroval);
    }

    std::vector<uint32_t> indices;
    std::vector<uint32_t> lengths;
    bool newstart = true;

    uint32_t sig_len = 0;
    uint32_t segment_len = 0;
    for (size_t i = 0; i < numints; i++) {
        if (inbuf[i] > zerothresh) {
            if (newstart) {
                segment_len = 0;
                indices.push_back(i);
                newstart = false;
            }
            segment_len++;
        } else if (!newstart) {
            lengths.push_back(segment_len);
            sig_len += segment_len;
            newstart = true;
        }
    }
    if (!newstart) {
        lengths.push_back(segment_len);
        sig_len += segment_len;
    }

    assert(indices.size() == lengths.size());
    dbg("%zu non-zero sections", indices.size());

    int32_t* sigs = (int32_t*)malloc(sig_len * sizeof(int32_t));
    size_t sig_needle = 0;
    for (size_t i = 0; i < indices.size(); i++) {
        for (size_t j = 0; j < lengths[i]; j++) {
            sigs[sig_needle] = inbuf[indices[i]+j];
            sig_needle++;
        }
    }

    delta_encode(sigs, sig_len);

    int32_t* percentiles = (int32_t*)malloc(sig_len * sizeof(int32_t));
    memcpy(percentiles, sigs, sig_len * sizeof(int32_t));

    // TODO: what if we want to look at less than the 99th percentile?
    std::nth_element(percentiles, percentiles + (sig_len * 99)/100, percentiles + sig_len);
    const int32_t p99 = percentiles[(sig_len*99)/100];
    std::nth_element(percentiles, percentiles + sig_len/100, percentiles + (sig_len * 99)/100);
    const int32_t p01 = percentiles[sig_len/100];
    dbg("Middle 98%% of deltas is [%d, %d]", p01, p99);

    free(percentiles);

    const uint16_t bound = std::max(std::abs(p99), std::abs(p01));
    // Ceiling of log2(bound), plus one for the sign bit
    const uint32_t bits_needed = (uint32_t)std::log2(bound) + 1 + 1;
    dbg("Bit-packing with %d bits", bits_needed);

    // marker is 00..0011...11 with bits_needed 1s
    const uint32_t ones = ~0;
    const uint32_t outlier_marker = ones >> (32 - bits_needed);

    uint32_t* zigzagged = (uint32_t*)malloc(sig_len * sizeof(uint32_t));
    // Don't reeeeally need 32 bits for these, just 16 + sign bit
    std::vector<int32_t> outliers;
    for (int i = 0; i < sig_len; i++) {
        if (std::abs(sigs[i]) > bound) {
            outliers.push_back(sigs[i]);
            zigzagged[i] = outlier_marker;
        } else {
            // zig-zag encode
            zigzagged[i] = (sigs[i] << 1) ^ (sigs[i] >> 31);
        }
    }
    free(sigs);
    dbg("%zu outliers encoded", outliers.size());

    uint32_t packed_bytes = simdpack_compressedbytes(sig_len, bits_needed);
    uint8_t* packed = (uint8_t*)malloc(packed_bytes);
    __m128i* endofbuf = simdpack_length(zigzagged, sig_len, (__m128i*)packed, bits_needed);
    packed_bytes = (endofbuf - (__m128i*)packed)*sizeof(__m128i);
    free(zigzagged);

    size_t offset = 0;

    struct header h;
    h.indexsize = indices.size();
    h.numints = numints;
    h.zeroval = zeroval;
    h.outliersize = outliers.size();
    h.siglen = sig_len;
    h.bitsneeded = bits_needed;

    offset = writeheader(outbuf, offset, h);

    writebig(outbuf, offset, indices.size()  * sizeof(uint32_t), indices.data());
    writebig(outbuf, offset, lengths.size()  * sizeof(uint32_t), lengths.data());
    writebig(outbuf, offset, outliers.size() * sizeof(int32_t),  outliers.data());
    writebig(outbuf, offset, packed_bytes, packed);

    free(packed);

    return offset;
}

void print_usage() {
    printf("Usage: compress [-b/i] in compressed\n");
    printf("  Compresses `in', a newline-separated list of numbers, into `compressed'");
    printf("  -b\tTakes in as a binary file of 16-bit integers\n");
    printf("  -i\tTakes in as a newline-separated list of text integers\n");
}

int main(const int argc, const char *argv[]) {
    size_t bytes;
    uint16_t* inbuf;
    std::string outfile;
    if (argv[1][0] == '-' && argv[1][1] == 'b' && argv[1][2] == '\0') {
        // binary mode
        if (argc != 4) {
            print_usage();
            return 1;
        }
        printf("Working in binary mode.\n");
        outfile = argv[3];
        printf("Compressing %s into %s.\n", argv[2], outfile.c_str());
        std::ifstream ifs(argv[2], std::ios::binary);
        std::vector<char> v = std::vector(std::istreambuf_iterator<char>(ifs), {});
        ifs.close();
        bytes = v.size();
        inbuf = (uint16_t*)v.data();
    } else if (argv[1][0] == '-' && argv[1][1] == 'i' && argv[1][2] == '\0') {
        // text int mode
        if (argc != 4) {
            print_usage();
            return 1;
        }

        printf("Working in text, integer mode.\n");
        outfile = argv[3];
        printf("Compressing %s into %s.\n", argv[2], outfile.c_str());
        std::ifstream ifs(argv[2]);
        std::vector<uint16_t> v(std::istream_iterator<uint16_t>(ifs), {});
        ifs.close();
        bytes = v.size();
        inbuf = v.data();
    } else {
        // text float moad
        if (argc != 3) {
            print_usage();
            return 1;
        }

        printf("Working in text, float mode.\n");
        outfile = argv[2];
        printf("Compressing %s into %s.\n", argv[1], outfile.c_str());
        std::ifstream ifs(argv[1]);
        std::vector<double> doubles(std::istream_iterator<double>(ifs), {});
        ifs.close();
        std::vector<uint16_t> v(doubles.begin(), doubles.end());
        bytes = v.size();
        inbuf = v.data();
    }

    char* outbuf = (char*)malloc(bytes * sizeof(char));
    std::ofstream ofs(outfile, std::ios::binary);
    if (!ofs) return 1; // fail

    auto start = std::chrono::steady_clock::now();
    int64_t outlen = deltacode(inbuf, bytes, outbuf);
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff = end - start;

    ofs.write(outbuf, outlen);
    ofs.close();
    free(outbuf);

    printf("Compressed in %f s. \U000026A1\n", diff.count());
    return 0;
}
