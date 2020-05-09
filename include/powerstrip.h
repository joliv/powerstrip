#ifndef POWERSTRIP_INCLUDE_COMPRESS_H_
#define POWERSTRIP_INCLUDE_COMPRESS_H_

#include <cstdint>
#include <cstddef>
#include <vector>

class Bitpacked {
  uint8_t* packed;
  uint32_t len;
  uint8_t bits;
  uint32_t bytes;

  Bitpacked(uint8_t* packed, uint32_t len, uint8_t bits, uint32_t bytes);

 public:
  [[nodiscard]] uint8_t GetBits() const;
  [[nodiscard]] uint32_t GetLen() const;
  static Bitpacked Pack(const std::vector<uint32_t> &xs, uint8_t bits);
  uint64_t Write(char* to) const;
  static Bitpacked Read(const char* from, uint64_t &offset);
  uint32_t* Unpack();
};

class Packed {
  Bitpacked signal;
  Bitpacked outliers;

  static int32_t* DeltaEncode(const uint16_t* xs, uint32_t len);
  static uint8_t BestBits(const int32_t* xs, uint32_t len);
  static uint32_t ZigZag(int32_t x);
  static int32_t UnZig(uint32_t x);
  Packed(const Bitpacked &signal, const Bitpacked &outliers);

 public:
  uint64_t Write(char* to) const;
  static Packed Read(const char* from, uint64_t &offset);
  static Packed Pack(const uint16_t* xs, uint32_t);
  uint16_t* Unpack(uint32_t* len);
};

class Stripped {
  uint32_t total_l;
  uint16_t* actives;
  uint32_t actives_l;
  uint32_t* indices;
  uint32_t* lengths;
  uint32_t segments;
  uint16_t floor;

  static int IsActive(uint16_t x, uint16_t floor);
 public:
  Stripped(uint32_t total_l,
           uint16_t* actives,
           uint32_t actives_l,
           uint32_t* indices,
           uint32_t* lengths,
           uint32_t segments,
           uint16_t floor);

 public:
  static Stripped Read(const char* from, uint64_t &offset, uint16_t* actives, uint32_t actives_l);
  static Stripped Strip(const uint16_t* xs, size_t len);
  uint64_t Unstrip(uint16_t* to);
  [[nodiscard]] uint16_t* GetActives() const;
  [[nodiscard]] uint32_t GetActivesL() const;
  uint64_t Write(char* to) const;
};

// Compress a block of len uint16s into out
uint64_t compress_block(const uint16_t* block, size_t len, char* out);
uint64_t decompress_block(const char* block, size_t len, uint16_t* out);

#endif //POWERSTRIP_INCLUDE_COMPRESS_H_
