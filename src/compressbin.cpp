#include <cstring>
#include <fstream>
#include <iostream>

#include "powerstrip.h"

void print_compression_stats(const uint64_t total_in,
                             const uint64_t total_out,
                             const std::chrono::time_point<std::chrono::steady_clock> start,
                             const std::chrono::time_point<std::chrono::steady_clock> end) {
  const double perc = ((double) total_out / (double) total_in) * 100;
  const std::chrono::duration<double> diff = end - start;
  std::cout << "Compressed to ";
  std::cout << std::fixed << std::setprecision(2) << perc;
  std::cout << "% of the original size";
  std::cout << " in " << diff.count() << " seconds \U000026A1" << std::endl;
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

  // ifstream is quite fast:
  // https://lemire.me/blog/2012/06/26/which-is-fastest-read-fread-ifstream-or-mmap/
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

  std::cout << "Compressing " << argv[1] << " into " << argv[2] << std::endl;

  char* i_buf = new char[INPUT_SIZE];
  char* o_buf = new char[OUTPUT_SIZE];

  uint64_t total_in = 0;
  uint64_t total_out = 0;
  const auto start = std::chrono::steady_clock::now();

  while (!i_f.eof()) {
    i_f.read(i_buf, INPUT_SIZE);
    const std::streamsize in_bytes = i_f.gcount();

    if (in_bytes % 2 != 0) {
      std::cout << "ERROR stream ends with half of a uint16" << std::endl;
      return -1;
    }

    total_in += in_bytes;
    const size_t len = in_bytes / sizeof(uint16_t);

    const int64_t out_bytes = compress_block((uint16_t*) i_buf, len, o_buf);
    total_out += out_bytes;

    o_f.write(reinterpret_cast<const char*>(&out_bytes), sizeof(out_bytes));
    o_f.write(o_buf, out_bytes);
  }

  delete[] i_buf;
  delete[] o_buf;

  const auto end = std::chrono::steady_clock::now();
  print_compression_stats(total_in, total_out, start, end);

  return 0;
}
