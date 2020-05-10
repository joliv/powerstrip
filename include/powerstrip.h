#ifndef POWERSTRIP_INCLUDE_POWERSTRIP_H_
#define POWERSTRIP_INCLUDE_POWERSTRIP_H_

#include <cstdint>

struct stripped {
  uint32_t segments;
  uint32_t* indices;
  uint32_t* lengths;
  uint32_t actives_l;
  uint16_t floor;
  uint32_t total_l;
};

struct bitpacked {
  uint8_t* packed;
  uint32_t len;
  uint8_t bits;
  uint32_t bytes;
};

struct packed {
  struct bitpacked signal;
  struct bitpacked outliers;
};

// Compress a block of len uint16s into out
uint64_t compress_block(const uint16_t* block, size_t len, char* out);
uint64_t decompress_block(const char* block, size_t len, uint16_t* out);

#endif //POWERSTRIP_INCLUDE_POWERSTRIP_H_
