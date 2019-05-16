#include <stdio.h>
#include <cstdint>
#include <fstream>
#include <vector>

#include "emmintrin.h" // for __m128, intel intrinsic
#include "simdcomp.h" // for bit packing

#include "decompress.hpp"

int64_t deltadecode(char* inbuf, size_t insize, uint16_t* outbuf) {
    // What follows could probably be a macro.
    size_t offset = 0;

    // Get structure lengths
    uint32_t* initiallen = reinterpret_cast<uint32_t*>(inbuf + offset);
    uint32_t indexlen = initiallen[0];
    offset += sizeof(uint32_t);

    uint32_t* totallen = reinterpret_cast<uint32_t*>(inbuf + offset);
    uint32_t total_len = totallen[0];
    offset += sizeof(uint32_t);

    uint16_t* zerovalbytes = reinterpret_cast<uint16_t*>(inbuf + offset);
    uint16_t zeroval = zerovalbytes[0];
    offset += sizeof(uint16_t);

    uint32_t* outlierlenbytes = reinterpret_cast<uint32_t*>(inbuf + offset);
    uint32_t outlier_len = outlierlenbytes[0];
    offset += sizeof(uint32_t);

    uint32_t* siglenbytes = reinterpret_cast<uint32_t*>(inbuf + offset);
    uint32_t sig_len = siglenbytes[0];
    offset += sizeof(uint32_t);

    uint8_t* bitneededbytes = reinterpret_cast<uint8_t*>(inbuf + offset);
    uint8_t bits_needed = bitneededbytes[0];
    offset += sizeof(uint8_t);

    uint32_t* indices = reinterpret_cast<uint32_t*>(inbuf + offset);
    offset += indexlen * sizeof(uint32_t);

    uint32_t* lengths = reinterpret_cast<uint32_t*>(inbuf + offset);
    offset += indexlen * sizeof(uint32_t);

    int32_t* outliers = reinterpret_cast<int32_t*>(inbuf + offset);
    offset += outlier_len * sizeof(int32_t);

    uint8_t* packed = reinterpret_cast<uint8_t*>(inbuf + offset);

    uint32_t* zigzagged = (uint32_t*)malloc(sig_len * sizeof(uint32_t));
    simdunpack_length((const __m128i*)packed, sig_len, zigzagged, bits_needed);

    // Might prefer to re-use memory here...decoding needs to be super-fast
    // int32_t* sigs = (int32_t*)malloc(sig_len * sizeof(int32_t));

    uint32_t outlier_needle = 0;
    const uint32_t ones = ~0;
    const uint32_t outlier_marker = ones >> (32 - bits_needed);
    uint32_t prev = 0;
    for (int i = 0; i < sig_len; i++) {
        // Is there a faster bit-twiddling check?
        if (zigzagged[i] == outlier_marker) {
            int32_t delta = outliers[outlier_needle];
            outlier_needle++;
            zigzagged[i] = prev + delta;
            prev = zigzagged[i];
        } else {
            int32_t delta = (zigzagged[i] >> 1) ^ -(zigzagged[i] & 1);
            zigzagged[i] = prev + delta;
            prev = zigzagged[i];
        }
    }
    // free(zigzagged);

    // Clear output
    for (uint32_t i = 0; i < total_len; i++) {
        outbuf[i] = zeroval;
    }

    size_t inneedle = 0;
    for (uint32_t i = 0; i < indexlen; i++) {
        for (uint32_t j = 0; j < lengths[i]; j++) {
            outbuf[indices[i]+j] = (uint16_t)zigzagged[inneedle];
            inneedle++; // A C++ guru would do this in one line
        }
    }
    free(zigzagged);

    return total_len * sizeof(uint16_t);
}

int main(const int argc, const char *argv[]) {
    std::ifstream ifs(argv[1], std::ios::binary);
    std::vector<char> inbuf(std::istreambuf_iterator<char>(ifs), {});
    ifs.close();

    // Allocate 20 times the space we got
    char* outbuf = (char*)malloc(sizeof(char) * 10000000);
    std::ofstream ofs(argv[2], std::ios::binary);
    if (!ifs || !ofs) return 1; // fail

    int64_t outlen = deltadecode(inbuf.data(), inbuf.size(), (uint16_t*)outbuf);
    ofs.write(outbuf, outlen);
    free(outbuf);
    ofs.close();

    printf("Printed to where it's s'possed to go! 🎉\n");
    return 0;
}
