#include "rapidcheck.h"

#include "powerstrip.h"

#include <cstring>
#include <vector>
#include <algorithm>

int main() {
  rc::check("round-trip yields the same numbers",
            [](const std::vector<uint16_t> &l0) {
              auto l1 = l0;
              char* comped = new char[OUTPUT_SIZE];
              auto* decomped = new uint16_t[l1.size()];
              const uint64_t len = compress_block(l1.data(), l1.size(), comped);
              const uint64_t size = decompress_block(comped, len, decomped);
              RC_ASSERT(l0.size() == size);
              for (size_t i = 0; i < size; i++) {
                RC_LOG() << l0[i] << "\t" << decomped[i] << std::endl;
              }
              RC_ASSERT(std::memcmp(l0.data(), decomped, size) == 0);

              delete[] comped;
              delete[] decomped;
            });

  return 0;
}