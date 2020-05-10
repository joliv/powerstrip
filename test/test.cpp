#include "rapidcheck.h"

#include "powerstrip.h"

#include <cstring>
#include <vector>
#include <algorithm>

int main() {
  rc::check("round-trip yields the same numbers",
            [](const std::vector<uint16_t> &l0) {
              auto l1 = l0;
              // is the 1000 byte headroom necessary?
              char* comped = new char[BLOCK_SIZE + sizeof(uint32_t)];
              auto* decomped = new uint16_t[l1.size()];
              const uint64_t len = compress_block(l1.data(), l1.size(), comped);
              std::cout << "len=" << len << std::endl;
              const uint64_t size = decompress_block(comped, len, decomped);
              std::cout << "size=" << len << std::endl;
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