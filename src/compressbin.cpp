#include <iostream> // for cin, cout
#include <cstring> // for memcpy

#include "../extern/zstd/lib/common/huf.h" // Huffman coding
#include "../extern/simdcomp/include/simdcomp.h" // for bit-packing

#include "common.h"
#include "powerstrip.h"

int main() {
  char* in_buf = new char[BLOCK_SIZE];
  char* out_buf = new char[BLOCK_SIZE];

  while (!std::cin.eof()) {
    std::cin.read(in_buf, BLOCK_SIZE);
    std::streamsize in_bytes = std::cin.gcount();

    if (in_bytes % 2 != 0) {
      std::cerr << "ERR stream ends with half of a uint16" << std::endl;
      return -1;
    }

    const size_t len = in_bytes / sizeof(uint16_t);

    const int64_t out_bytes = compress_block((uint16_t*) in_buf, len, out_buf);

    std::cout.write(reinterpret_cast<const char*>(&out_bytes), sizeof(out_bytes));
    std::cout.write(out_buf, out_bytes);
  }

  delete[] in_buf;
  delete[] out_buf;

  std::cerr << "Hello, World!" << std::endl;
  std::cerr << HUF_isError(5) << std::endl;
  std::cerr << simdpack_compressedbytes(5, 5) << std::endl;

  return 0;
}
