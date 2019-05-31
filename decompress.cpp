#include <stdio.h>
#include <cstdint>
#include <fstream>
#include <vector>

#include "emmintrin.h" // for __m128, intel intrinsic
#include "simdcomp.h" // for bit packing

#include "decompress.hpp"

#ifdef DEBUG
#define dbg(...) do {\
    printf("  ");\
    printf(__VA_ARGS__);\
    printf("\n");\
} while(0)
#else
#define dbg(...) do {} while(0);
#endif

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

size_t get_numints(char* inbuf) {
    struct header h;
    readheader(inbuf, &h);

    return h.numints;
}

int64_t deltadecode(char* inbuf, size_t insize, uint16_t* outbuf) {
    // What follows could probably be a macro.
    size_t offset = 0;

    struct header h;

    offset += readheader(inbuf, &h);
    dbg("%d segments, %d outliers", h.indexsize, h.outliersize);

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

void print_usage() {
    printf("Usage: decompress [-b/i] in decompressed\n");
    printf("  Decompresses `in' into `decompressed', a newline-separated list of integers");
    printf("  -b\tWrites a binary file of 16-bit integers\n");
    printf("  -i\tPerforms the default action. This flag has no effect.\n");
}

enum Mode { mode_binary, mode_ints, mode_floats, mode_default };

int main(const int argc, const char *argv[]) {
    Mode mode;
    if (argv[1][0] != '-') {
        mode = mode_default;
        printf("Working in text mode.");
    } else if (argv[1][1] == 'b' && argv[1][2] == '\0') {
        mode = mode_binary;
        printf("Working in binary mode.");
    } else if (argv[1][1] == 'i' && argv[1][2] == '\0') {
        mode = mode_ints;
        printf("Working in text mode.");
    } else {
        printf("Unknown flag.\n");
        print_usage();
        exit(1);
    }
    std::ifstream ifs(mode == mode_default ? argv[1] : argv[2], std::ios::binary);
    std::vector<char> inbuf(std::istreambuf_iterator<char>(ifs), {});
    ifs.close();

    const size_t numints = get_numints(inbuf.data());
    uint16_t* outbuf = (uint16_t*)malloc(sizeof(uint16_t) * numints);
    dbg("Decompressing %zu ints", numints);
    std::ofstream ofs(mode == mode_default ? argv[2] : argv[3], std::ios::binary);
    if (!ifs || !ofs) return 1; // fail

    auto start = std::chrono::steady_clock::now();
    int64_t outlen = deltadecode(inbuf.data(), inbuf.size(), (uint16_t*)outbuf);
    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> diff = end - start;

    if (mode == mode_binary) {
        // binary mode
        ofs.write((char*)outbuf, outlen);
    } else {
        // text int mode
        for (int i = 0; i < outlen/sizeof(uint16_t); i++) {
            ofs << outbuf[i] << '\n';
        }
    }
    free(outbuf);
    ofs.close();

    printf("Decompressed in %f s. \U000026A1\n", diff.count());
    return 0;
}
