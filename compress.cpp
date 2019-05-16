#include <stdio.h>
#include <cmath> // for abs, log2
#include <cstdint> // for fixed ints like uint16_t
#include <fstream>
#include <utility>
#include <vector> // for vector
#include <cstring>
#include <algorithm> // for max

#include "emmintrin.h" // for __m128, intel intrinsic
#include "vendor/simdcomp/include/simdcomp.h" // for bit packing

#include "compress.hpp"

#define writesmall(buf, offset, type, x) do {\
    type* interp = reinterpret_cast<type*>(buf+offset);\
    interp[0] = x;\
    offset += sizeof(type);\
} while(0)

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

int64_t deltacode(uint16_t* inbuf, size_t insize, char* outbuf) {
    const size_t numints = insize / 2; // 2 bytes per uint16

    // const uint16_t zeroval = 0, zerothresh = 0;
    uint16_t zeroval = get_zeroval(inbuf, numints);
    uint16_t zerothresh = zeroval + 5; // Sure, let's say +5W for now

    // No floor exists.
    if (zeroval == UINT16_MAX) {
        zeroval = 0;
        zerothresh = 0;
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
    // printf("\n99th percentile is %d\n", percentiles[(sig_len*99)/100]);
    std::nth_element(percentiles, percentiles + sig_len/100, percentiles + (sig_len * 99)/100);
    const int32_t p01 = percentiles[sig_len/100];
    // printf("\n1st percentile is %d\n", percentiles[sig_len/100]);

    free(percentiles);

    const uint16_t bound = std::max(std::abs(p99), std::abs(p01));
    // Ceiling of log2(bound), plus one for the sign bit
    const uint32_t bits_needed = (uint32_t)std::log2(bound) + 1 + 1;
    // printf("\nWe'll need %d bits for this.\n", bits_needed);

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

    uint32_t packed_bytes = simdpack_compressedbytes(sig_len, bits_needed);
    uint8_t* packed = (uint8_t*)malloc(packed_bytes);
    __m128i* endofbuf = simdpack_length(zigzagged, sig_len, (__m128i*)packed, bits_needed);
    packed_bytes = (endofbuf - (__m128i*)packed)*sizeof(__m128i);
    free(zigzagged);

    size_t offset = 0;

    writesmall(outbuf, offset, uint32_t, indices.size());
    writesmall(outbuf, offset, uint32_t, numints);
    writesmall(outbuf, offset, uint16_t, zeroval);
    writesmall(outbuf, offset, uint32_t, outliers.size());
    writesmall(outbuf, offset, uint32_t, sig_len);
    writesmall(outbuf, offset, uint8_t,  bits_needed);

    writebig(outbuf, offset, indices.size()  * sizeof(uint32_t), indices.data());
    writebig(outbuf, offset, lengths.size()  * sizeof(uint32_t), lengths.data());
    writebig(outbuf, offset, outliers.size() * sizeof(int32_t),  outliers.data());
    writebig(outbuf, offset, packed_bytes, packed);

    free(packed);

    return offset;
}

int main(const int argc, const char *argv[]) {
    std::ifstream ifs(argv[1], std::ios::binary);
    std::vector<char> inbuf(std::istreambuf_iterator<char>(ifs), {});
    ifs.close();

    char* outbuf = (char*)malloc(inbuf.size() * sizeof(char));
    std::ofstream ofs(argv[2], std::ios::binary);
    if (!ifs || !ofs) return 1; // fail

    int64_t outlen = deltacode((uint16_t*)inbuf.data(), inbuf.size(), outbuf);
    ofs.write(outbuf, outlen);
    ofs.close();
    free(outbuf);

    printf("Printed to where it's s'possed to go! 🎉\n");
    return 0;
}
