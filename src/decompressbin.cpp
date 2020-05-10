#include <iostream>
#include <fstream>

#include "common.h"
#include "powerstrip.h"

int main(const int argc, const char** argv) {
  // TODO replace with help text
  if (argc != 3) {
    std::cerr << "ERR insufficient args" << std::endl;
    return -1;
  }

  // ifstream is pretty fast!
  // https://lemire.me/blog/2012/06/26/which-is-fastest-read-fread-ifstream-or-mmap/
  std::ifstream i_f(argv[1], std::ios::binary);
  std::ofstream o_f(argv[2], std::ios::binary);

  auto* in_buf = new char[BLOCK_SIZE];
  auto* out_buf = new uint16_t[BLOCK_SIZE / sizeof(uint16_t)];

  while (!i_f.eof()) {
    uint64_t block_len;
    i_f.read(reinterpret_cast<char*>(&block_len), sizeof(block_len));
    if (i_f.gcount() == 0) break;

    i_f.read(in_buf, block_len);

    const uint64_t out_len = decompress_block(in_buf, block_len, out_buf);
    o_f.write(reinterpret_cast<const char*>(out_buf), out_len * sizeof(uint16_t));
  }

  delete[] in_buf;
  delete[] out_buf;

  return 0;
}