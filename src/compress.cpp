#include <iostream> // for cin, cout

#include "../extern/zstd/lib/common/huf.h" // Huffman coding
#include "../extern/simdcomp/include/simdcomp.h" // for bit-packing

// 128 KB = HUF_BLOCKSIZE_MAX
// TODO allow larger sizes but if it's too much we just fall back
#define BLOCK_SIZE (128 * 1024 * 4)

uint64_t compress_block(uint16_t* block, size_t len, char* out) {
  return 0;
}

int main() {
  char* in_buf = new char[BLOCK_SIZE];
  char* out_buf = new char[BLOCK_SIZE];

  while (!std::cin.eof()) {
    std::cin.read(in_buf, BLOCK_SIZE);
    std::streamsize in_len = std::cin.gcount();

    const int64_t out_len = compress_block((uint16_t*) in_buf, in_len, out_buf);

    std::cout.write(reinterpret_cast<const char*>(&out_len), sizeof(out_len));
    std::cout.write(out_buf, out_len);
  }

  delete[] in_buf;
  delete[] out_buf;

  std::cout << "Hello, World!" << std::endl;
  std::cout << HUF_isError(5) << std::endl;
  std::cout << simdpack_compressedbytes(5, 5) << std::endl;

  return 0;
}
