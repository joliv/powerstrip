#include <iostream> // for cin, cout
#include <cstring> // for memcpy
#include <fstream> // for if,ofstream

#include "../extern/zstd/lib/common/huf.h" // Huffman coding
#include "../extern/simdcomp/include/simdcomp.h" // for bit-packing

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

  std::cerr << "Compressing from " << argv[1] << " into " << argv[2] << std::endl;

  char* i_buf = new char[BLOCK_SIZE];
  char* o_buf = new char[BLOCK_SIZE + sizeof(uint32_t)]; // One number of headroom for Huffman length

  while (!i_f.eof()) {
    i_f.read(i_buf, BLOCK_SIZE);
    std::streamsize in_bytes = i_f.gcount();

    if (in_bytes % 2 != 0) {
      std::cerr << "ERR stream ends with half of a uint16" << std::endl;
      return -1;
    }

    const size_t len = in_bytes / sizeof(uint16_t);

    const int64_t out_bytes = compress_block((uint16_t*) i_buf, len, o_buf);

    o_f.write(reinterpret_cast<const char*>(&out_bytes), sizeof(out_bytes));
    o_f.write(o_buf, out_bytes);
  }

  delete[] i_buf;
  delete[] o_buf;

  std::cerr << "Hello, World!" << std::endl;

  return 0;
}
