#ifndef POWERSTRIP_INCLUDE_POWERSTRIP_H_
#define POWERSTRIP_INCLUDE_POWERSTRIP_H_

#include <cstdint>

// Based on experimental results.
// One thing is: once power-stripping is done, if we're above HUF_BLOCKSIZE_MAX, we won't Huffman-code.
#define INPUT_SIZE (128 * 1024 * 4)
#define OUTPUT_SIZE (INPUT_SIZE + sizeof(uint32_t))

// Compress a block of 'len' uint16s into 'out'
// ensure:
//   len * sizeof(uint16_t) <= INPUT_SIZE
//   out is an allocated block of at least OUTPUT_SIZE bytes
// returns: the size of the compressed data, in bytes
uint64_t compress_block(const uint16_t* block, size_t len, char* out);

// Decompress a block of 'len' bytes of Powerstripped data into 'out'
// ensure:
//   out is an allocated block of at least INPUT_SIZE bytes
// returns: the size of the uncompressed data, in bytes
uint64_t decompress_block(const char* block, size_t len, uint16_t* out);

#endif //POWERSTRIP_INCLUDE_POWERSTRIP_H_
