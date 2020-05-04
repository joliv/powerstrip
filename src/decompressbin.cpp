#include <iostream>

#include "common.h"
#include "powerstrip.h"

int main() {
  auto* in_buf = new char[BLOCK_SIZE];
  auto* out_buf = new uint16_t[BLOCK_SIZE / sizeof(uint16_t)];

  while (!std::cin.eof()) {
    uint64_t block_len;
    std::cin.read(reinterpret_cast<char*>(&block_len), sizeof(block_len));
    if (std::cin.gcount() == 0) break;

    std::cin.read(in_buf, block_len);

    const uint64_t out_len = decompress_block(in_buf, block_len, out_buf);
    std::cout.write(reinterpret_cast<const char*>(out_buf), out_len * sizeof(uint16_t));
  }

  delete[] in_buf;
  delete[] out_buf;

  return 0;
}