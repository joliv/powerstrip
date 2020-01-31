#include <stdio.h>   // for printf
#include <cmath>     // for abs, log2
#include <cstdint>   // for fixed ints like uint16_t
#include <fstream>   // for ifstream
#include <utility>   // I don't remember
#include <vector>    // for vector
#include <algorithm> // for max
#include <string>    // for std::string
#include <chrono>    // timing, for debugging

#include "emmintrin.h" // for __m128, intel intrinsic
#include "vendor/simdcomp/include/simdcomp.h" // for bit packing

#include "vendor/zstd/lib/common/huf.h" // Huffman coding

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
    memcpy(buf + offset, x, size);\
    offset += size;\
} while(0)

#define OUTLIER_BITS 17

// One pass to find the noise floor, basically just the mode
uint16_t get_zeroval(const uint16_t* inbuf, const size_t numints) {
    uint16_t common = 0;
    size_t common_count = 0;
    size_t* histogram = (size_t*)calloc(UINT16_MAX+1, sizeof(size_t));
    for (size_t i = 0; i < numints; i++) {
        histogram[inbuf[i]]++;
        if (histogram[inbuf[i]] > common_count) {
            common = inbuf[i];
            common_count = histogram[inbuf[i]];
        }
    }
    free(histogram);
    dbg("Floor of %d is %zu/%zu total", common, common_count, numints);
    // if (10 * common_count > numints) {
        // return common;
    // }
    return UINT16_MAX; // no floor with >10% of samples found
}

double measure_err(const uint16_t* inbuf, const size_t numints, const uint16_t zeroval, const uint16_t zerothresh) {
    double sumerr = 0;
    for (size_t i = 0; i < numints; i++) {
        if (inbuf[i] <= zerothresh) {
            sumerr += std::abs(inbuf[i] - zeroval);
        }
    }
    return sumerr / (double)numints;
}

// Simple delta encoding. See how lemire/simdcomp does this fast?
// (Nevermind: the answer is SIMD, of course)
void delta_encode(int32_t* xs, uint32_t len) {
    int32_t prev = 0;
    for (uint32_t i = 0; i < len; i++) {
        int32_t tmpprev = xs[i];
        xs[i] -= prev;
        prev = tmpprev;
    }
}

uint32_t best_bits(int32_t* xs, uint32_t len) {
    // Find the best number of bits to use for non-outliers
    size_t* bitcounts = (size_t*)calloc(17, sizeof(size_t));

    for (size_t i = 0; i < len; i++) {
        // assert((size_t)std::log2(std::abs(xs[i])) + 1 + 1 <= 16);
        if (xs[i] == 0) { // protect against log2(0) = -inf
            bitcounts[1]++;
        } else {
            // Ceiling of log2(bound), plus one for the sign bit
            bitcounts[(size_t)std::log2(std::abs(xs[i])) + 1 + 1]++;
        }
    }

    size_t bestbits = 0;
    size_t bestbittotal = SIZE_MAX;
    dbg("bits per int/total bits");
    for (size_t i = 1; i <= 16; i++) {
        bitcounts[i] += bitcounts[i-1];
        // 17 bytes for the outliers because 16 + 1 sign but
        size_t totalbits = i * bitcounts[i] + OUTLIER_BITS * (len - bitcounts[i]);
        dbg("  %zu with %zu bits: %zu packed + %zu outliers", i, totalbits, i * bitcounts[i], OUTLIER_BITS * (len - bitcounts[i]));
        if (totalbits < bestbittotal) {
            bestbittotal = totalbits;
            bestbits = i;
        }
    }

    free(bitcounts);

    // Ceiling of log2(bound), plus one for the sign bit
    return (uint32_t)bestbits;
}

struct header {
    uint32_t indexsize;
    uint32_t numints;
    uint16_t zeroval;
    uint32_t outliersize;
    uint32_t outlierlen;
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
    writesmall(outbuf, offset, uint32_t, h.outlierlen);
    writesmall(outbuf, offset, uint32_t, h.siglen);
    writesmall(outbuf, offset, uint8_t,  h.bitsneeded);
    dbg("  Header: %zu",
      sizeof(uint32_t) + sizeof(uint32_t) + sizeof(uint16_t)
      + sizeof(uint32_t) + sizeof(uint32_t) + sizeof(uint32_t)
      + sizeof(uint8_t));
    return offset;
}

uint32_t bitpack(uint32_t* in, uint8_t** out, uint32_t len, uint32_t bits) {
    uint32_t packed_bytes = simdpack_compressedbytes(len, bits);
    uint8_t* packed = (uint8_t*)malloc(packed_bytes);
    __m128i* endofbuf = simdpack_length(in, len, (__m128i*)packed, bits);
    packed_bytes = (endofbuf - (__m128i*)packed)*sizeof(__m128i);
    *out = packed;
    return packed_bytes;
}

static inline uint32_t zigzag(int32_t x) {
    return (x << 1) ^ (x >> 31);
}

int64_t deltacode(uint16_t* inbuf, size_t insize, char* outbuf) {
    const size_t numints = insize / 2; // 2 bytes per uint16
    dbg("Encoding %zu ints", numints);

    uint16_t zeroval = get_zeroval(inbuf, numints);
    uint16_t zerothresh = zeroval + 5; // Sure, let's say +5W for now

    // No floor exists.
    if (zeroval == UINT16_MAX) {
        zeroval = 0;
        zerothresh = 0;
        dbg("No floor found. Falling back to lossless");
        dbg("RMSE will be 0 watts");
    } else {
        dbg("Floor found at %d watts", zeroval);
        dbg("RMSE will be %f watts", measure_err(inbuf, numints, zeroval, zerothresh));
    }

    std::vector<uint32_t> indices;
    std::vector<uint32_t> lengths;
    bool newstart = true;

    // Get the locations of the significant sections
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
    dbg("%zu non-zero section(s) with sum length %d", indices.size(), sig_len);

    // Pack the significant sections into `sigs`
    int32_t* sigs = (int32_t*)malloc(sig_len * sizeof(int32_t));
    size_t sig_needle = 0;
    for (size_t i = 0; i < indices.size(); i++) {
        for (size_t j = 0; j < lengths[i]; j++) {
            sigs[sig_needle] = inbuf[indices[i]+j];
            sig_needle++;
        }
    }

    // Delta-encode `sigs` in-place
    delta_encode(sigs, sig_len);

    const uint32_t bits_needed = best_bits(sigs, sig_len);

    const uint32_t ones = ~0;
    // The marker is 00..0011...11 with bits_needed 1s
    // Would probably work with 11...11 too but I'm scared to try
    const uint32_t outlier_marker = ones >> (32 - bits_needed);
    // Everything below that top spot is fair game
    const uint32_t bound = outlier_marker - 1;
    dbg("Bit-packing with %d bits. All > %u are outliers", bits_needed, bound);

    uint32_t* zigzagged = (uint32_t*)malloc(sig_len * sizeof(uint32_t));
    std::vector<uint32_t> outliers;
    for (int i = 0; i < sig_len; i++) {
        uint32_t zigged = zigzag(sigs[i]);
        if (zigged > bound) {
            outliers.push_back(zigged);
            zigzagged[i] = outlier_marker;
        } else {
            // zig-zag encode
            zigzagged[i] = zigged;
        }
    }
    free(sigs);
    dbg("%zu outliers encoded", outliers.size());

    // Pack outliers into 17 bits
    uint8_t* packed_outliers;
    uint32_t outlier_bytes = bitpack(outliers.data(), &packed_outliers, outliers.size(), OUTLIER_BITS);

    // Pack zigzagged into bits_needed bits
    uint8_t* packed_sigs;
    uint32_t packed_bytes = bitpack(zigzagged, &packed_sigs, sig_len, bits_needed);
    free(zigzagged);

    size_t offset = 0;

    struct header h;
    h.indexsize = indices.size();
    h.numints = numints;
    h.zeroval = zeroval;
    h.outliersize = outlier_bytes;
    h.outlierlen = outliers.size();
    h.siglen = sig_len;
    h.bitsneeded = bits_needed;

    size_t output_size = sizeof(uint32_t) + sizeof(uint32_t) + sizeof(uint16_t)
      + sizeof(uint32_t) + sizeof(uint32_t) + sizeof(uint32_t)
      + sizeof(uint8_t)
      + indices.size() * sizeof(uint32_t)
      + lengths.size() * sizeof(uint32_t)
      + outlier_bytes
      + packed_bytes;

    char* before_compression = (char*)malloc(output_size);

    dbg("Output sizes, in bytes");

    offset = writeheader(before_compression, offset, h);

    dbg("  Indices: %zu", indices.size() * sizeof(uint32_t));
    writebig(before_compression, offset, indices.size() * sizeof(uint32_t), indices.data());
    dbg("  Lengths: %zu", lengths.size() * sizeof(uint32_t));
    writebig(before_compression, offset, lengths.size() * sizeof(uint32_t), lengths.data());
    dbg("  Outliers: %u", outlier_bytes);
    writebig(before_compression, offset, outlier_bytes, packed_outliers);
    dbg("  Packed data: %u", packed_bytes);
    writebig(before_compression, offset, packed_bytes, packed_sigs);
    dbg("  Total before Huffman: %zu", offset);

    free(packed_sigs);

    assert(offset == output_size);

    size_t x = HUF_compress(outbuf + sizeof(uint32_t), insize, before_compression, output_size);
    if (HUF_isError(x)) {
        x = 0; // let's just not compress
        printf("OH NO, ERRROROROROROROR\n");
    } else {
        printf("no error.\n");
    }
    dbg("  Total after Huffman: %zu", x);
    size_t zero = 0;
    if (x == 0) {
        // the data is uncompressible
        memcpy(outbuf + sizeof(uint32_t), before_compression, output_size);
        writesmall(outbuf, zero, uint32_t, 0);
    } else {
        writesmall(outbuf, zero, uint32_t, output_size);
    }

    return x;
}

void print_usage() {
    printf("Usage: compress [-b/i] in compressed\n");
    printf("  Compresses `in', a newline-separated list of numbers, into `compressed'\n");
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
        inbuf = (uint16_t*)malloc(bytes);
        memcpy(inbuf, v.data(), bytes);
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
        bytes = v.size() * sizeof(uint16_t);
        inbuf = (uint16_t*)malloc(bytes);
        memcpy(inbuf, v.data(), bytes);
    } else {
        // text float mode
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
        bytes = v.size() * sizeof(uint16_t);
        inbuf = (uint16_t*)malloc(bytes);
        memcpy(inbuf, v.data(), bytes);
    }

    // 100 bytes for a little headroom. This guess could probably be improved.
    char* outbuf = (char*)malloc(bytes * sizeof(char) + 100);
    std::ofstream ofs(outfile, std::ios::binary);
    if (!ofs) return 1; // fail

    auto start = std::chrono::steady_clock::now();
    int64_t outlen = deltacode(inbuf, bytes, outbuf);
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff = end - start;

    ofs.write(outbuf, outlen);
    ofs.close();
    free(inbuf);
    free(outbuf);

    printf("Compressed in %f s. \U000026A1\n", diff.count());
    return 0;
}
