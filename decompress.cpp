#include <stdio.h>
#include <cstdint>
#include <fstream>
#include <vector>

#include "emmintrin.h" // for __m128, intel intrinsic
#include "simdcomp.h" // for bit packing

#include "decompress.hpp"

struct header {
    uint32_t indexsize;
    uint32_t numints;
    uint16_t zeroval;
    uint32_t outliersize;
    uint32_t siglen;
    uint8_t  bitsneeded;
};

#define readsmall(buf, offset, type, s, entry) do {\
    type* x = reinterpret_cast<type*>(buf + offset);\
    s->entry = x[0];\
    offset += sizeof(type);\
} while(0)

size_t readheader(char* inbuf, struct header* h) {
    size_t offset = 0;

    readsmall(inbuf, offset, uint32_t, h, indexsize);
    readsmall(inbuf, offset, uint32_t, h, numints);
    readsmall(inbuf, offset, uint16_t, h, zeroval);
    readsmall(inbuf, offset, uint32_t, h, outliersize);
    readsmall(inbuf, offset, uint32_t, h, siglen);
    readsmall(inbuf, offset, uint8_t,  h, bitsneeded);

    return offset;
}

int64_t deltadecode(char* inbuf, size_t insize, uint16_t* outbuf) {
    // What follows could probably be a macro.
    size_t offset = 0;

    struct header h;

    offset += readheader(inbuf, &h);

    uint32_t* indices = reinterpret_cast<uint32_t*>(inbuf + offset);
    offset += h.indexsize * sizeof(uint32_t);

    uint32_t* lengths = reinterpret_cast<uint32_t*>(inbuf + offset);
    offset += h.indexsize * sizeof(uint32_t);

    int32_t* outliers = reinterpret_cast<int32_t*>(inbuf + offset);
    offset += h.outliersize * sizeof(int32_t);

    uint8_t* packed = reinterpret_cast<uint8_t*>(inbuf + offset);

    uint32_t* zigzagged = (uint32_t*)malloc(h.siglen * sizeof(uint32_t));
    simdunpack_length((const __m128i*)packed, h.siglen, zigzagged, h.bitsneeded);

    uint32_t outlier_needle = 0;
    const uint32_t ones = ~0;
    const uint32_t outlier_marker = ones >> (32 - h.bitsneeded);
    uint32_t prev = 0;
    for (int i = 0; i < h.siglen; i++) {
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

    // Clear output
    for (uint32_t i = 0; i < h.numints; i++) {
        outbuf[i] = h.zeroval;
    }

    size_t inneedle = 0;
    for (uint32_t i = 0; i < h.indexsize; i++) {
        for (uint32_t j = 0; j < lengths[i]; j++) {
            outbuf[indices[i]+j] = (uint16_t)zigzagged[inneedle];
            inneedle++; // A C++ guru would do this in one line
        }
    }
    free(zigzagged);

    return h.numints * sizeof(uint16_t);
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

    printf("Printed to where it's s'possed to go! \U000026A1\n");
    return 0;
}
