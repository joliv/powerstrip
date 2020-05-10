#ifndef POWERSTRIP_INCLUDE_POWERSTRIP_H_
#define POWERSTRIP_INCLUDE_POWERSTRIP_H_

#include <cstdint>

// 128 KB = HUF_BLOCKSIZE_MAX
// TODO I'm not sure why we multiply by 4?
#define BLOCK_SIZE (128 * 1024 * 4)

// Compress a block of len uint16s into out
uint64_t compress_block(const uint16_t* block, size_t len, char* out);
uint64_t decompress_block(const char* block, size_t len, uint16_t* out);

#endif //POWERSTRIP_INCLUDE_POWERSTRIP_H_
