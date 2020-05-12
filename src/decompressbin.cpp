#include <fstream>
#include <iostream>

#include "powerstrip.h"

void print_compression_stats(const std::chrono::time_point<std::chrono::steady_clock> start,
                             const std::chrono::time_point<std::chrono::steady_clock> end) {
  const std::chrono::duration<double> diff = end - start;
  std::cout << "Decompressed in " << std::fixed << std::setprecision(2) << diff.count() << " seconds \U000026A1"
            << std::endl;
}

void print_usage(const char* exec) {
  std::cout << "Usage:" << std::endl;
  std::cout << "  " << exec << " input.bin output.pstrip" << std::endl;
}

int main(const int argc, const char** argv) {
  if (argc != 3) {
    std::cout << "ERROR insufficient args" << std::endl;
    print_usage(argv[0]);
    return -1;
  }

  std::ifstream i_f(argv[1], std::ios::binary);
  if (!i_f) {
    std::cout << "Failed to open " << argv[1] << std::endl;
    return -1;
  }

  std::ofstream o_f(argv[2], std::ios::binary | std::ios::trunc);
  if (!o_f) {
    std::cout << "Failed to open " << argv[2] << std::endl;
    return -1;
  }

  std::cout << "Decompressing " << argv[1] << " into " << argv[2] << std::endl;

  auto* in_buf = new char[OUTPUT_SIZE];
  auto* out_buf = new uint16_t[INPUT_SIZE / sizeof(uint16_t)];

  const auto start = std::chrono::steady_clock::now();

  while (!i_f.eof()) {
    uint64_t block_len;
    i_f.read(reinterpret_cast<char*>(&block_len), sizeof(block_len));
    if (i_f.gcount() == 0) break;

    i_f.read(in_buf, block_len);

    const uint64_t out_len = decompress_block(in_buf, block_len, out_buf);
    o_f.write(reinterpret_cast<const char*>(out_buf), out_len * sizeof(uint16_t));
  }

  const auto end = std::chrono::steady_clock::now();
  print_compression_stats(start, end);

  delete[] in_buf;
  delete[] out_buf;

  return 0;
}