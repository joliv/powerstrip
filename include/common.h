#ifndef POWERSTRIP_INCLUDE_COMMON_H_
#define POWERSTRIP_INCLUDE_COMMON_H_

// 128 KB = HUF_BLOCKSIZE_MAX
// TODO allow larger sizes but if it's too much we just fall back
#define BLOCK_SIZE (128 * 1024 * 4)

// TODO why do we still have a "common" header?
#ifndef NDEBUG // cross-plat!
#define dbg(...) do {\
    printf("  ");\
    printf(__VA_ARGS__);\
    printf("\n");\
} while(0)
#else
#define dbg(...) do {} while(0);
#endif

#endif //POWERSTRIP_INCLUDE_COMMON_H_
