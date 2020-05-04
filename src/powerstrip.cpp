#include <cstring> // for memcpy
#include <utility>
#include <vector>
#include <cinttypes> // for printf macros
#include <powerstrip.h>
#include <numeric>
#include <array>
#include <limits>
#include <common.h>

#include "../extern/zstd/lib/common/huf.h" // Huffman coding
#include "../extern/simdcomp/include/simdcomp.h" // for bit-packing

#define WINDOW 0

#define OUTLIER_BITS 17

// TODO namespace?
// TODO use arena/workspace?
// TODO replace mallocs with new?

#define write_sm(buf, offset, type, x) do {\
  (reinterpret_cast<type*>(buf+offset))[0] = x;\
  offset += sizeof(type);\
} while(0)

#define write_lg(buf, offset, size, x) do {\
  memcpy(buf + offset, x, size);\
  offset += size;\
} while(0)

#define read_sm(buf, offset, type, x) do {\
  x = (reinterpret_cast<const type*>(buf + offset))[0];\
  offset += sizeof(type);\
} while(0)

uint64_t write(const Stripped &s, const Packed &p, char* to) {
  size_t offset = 0;

  offset += p.Write(to + offset);
  offset += s.Write(to + offset);

  return offset;
}

uint64_t compress_block(const uint16_t* block, const size_t len, char* out) {
  dbg("~");
  dbg("comp:  len=%zu", len);
  auto s = Stripped::Strip(block, len);
  auto p = Packed::Pack(s.GetActives(), s.GetActivesL());

  return write(s, p, out);
}

uint64_t decompress_block(const char* block, const size_t len, uint16_t* out) {
  uint64_t offset = 0;
  Packed p = Packed::Read(block, offset);
  uint32_t unpacked_l;
  uint16_t* unpacked = p.Unpack(&unpacked_l);
  Stripped s = Stripped::Read(block, offset, unpacked, unpacked_l);
  return s.Unstrip(out);
}

int Stripped::IsActive(const uint16_t x, const uint16_t floor) {
  return x > floor + WINDOW || x < floor - WINDOW;
}

Stripped Stripped::Strip(const uint16_t* xs, const size_t len) {
  uint16_t floor = 0; // TODO actually find this

  bool was_active = false;
  uint32_t segment_l = 0;
  uint32_t segments = 0;
  uint32_t actives_l = 0;

  if (len >= 1) dbg("strip: xs[0]=%" PRIu16 " len=%zu", xs[0], len);

  auto* indices = static_cast<uint32_t*>(std::malloc(len * sizeof(uint32_t)));
  auto* lengths = static_cast<uint32_t*>(std::malloc(len * sizeof(uint32_t)));
  auto* actives = static_cast<uint16_t*>(std::malloc(len * sizeof(uint16_t)));

  for (size_t i = 0; i < len; i++) {
    // TODO consider inverting the nesting so it's more like an FSM
    if (IsActive(xs[i], floor)) {
      if (!was_active) {
        segment_l = 0;
        indices[segments++] = i;
        was_active = true;
      }
      actives[actives_l++] = xs[i];
      segment_l++;
    } else if (was_active) {
      lengths[segments - 1] = segment_l;
      was_active = false;
    }
  }

  if (was_active) {
    lengths[segments - 1] = segment_l;
  }

  dbg("strip: actives_l=%" PRIu32 " segments=%" PRIu32, actives_l, segments);

  return Stripped(len, actives, actives_l, indices, lengths, segments, floor);
}

uint64_t Stripped::Unstrip(uint16_t* to) {
  for (uint32_t i = 0; i < total_l; i++) {
    to[i] = floor;
  }

  uint32_t needle = 0;
  for (uint32_t i = 0; i < segments; i++) {
    for (uint32_t j = 0; j < lengths[i]; j++) {
      dbg("ustrp: to[indices[%" PRIu32 "]+%" PRIu32 "]=actives[%" PRIu32 "]=%" PRIu16, i, j, needle, actives[needle]);
      to[indices[i] + j] = actives[needle++];
    }
  }

  return total_l;
}

Stripped::Stripped(uint32_t total_l,
                   uint16_t* actives,
                   uint32_t actives_l,
                   uint32_t* indices,
                   uint32_t* lengths,
                   uint32_t segments,
                   uint16_t floor)
    : total_l(total_l),
      actives(actives),
      actives_l(actives_l),
      indices(indices),
      lengths(lengths),
      segments(segments),
      floor(floor) {}

uint16_t* Stripped::GetActives() const {
  return actives;
}

// TODO hmm a little ugly
uint32_t Stripped::GetActivesL() const {
  return actives_l;
}

uint64_t Stripped::Write(char* to) const {
  uint64_t offset = 0;

  dbg("write: total_l=%" PRIu32 " floor=%" PRIu16 " segments=%" PRIu32, total_l, floor, segments);
  write_sm(to, offset, uint32_t, total_l);
  write_sm(to, offset, uint16_t, floor);
  write_sm(to, offset, uint32_t, segments);
  write_lg(to, offset, segments * sizeof(indices[0]), indices);
  write_lg(to, offset, segments * sizeof(lengths[0]), lengths);

  return offset;
}

Stripped Stripped::Read(const char* from, uint64_t &offset, uint16_t* actives, uint32_t actives_l) {
  uint32_t total_l;
  read_sm(from, offset, uint32_t, total_l);

  uint16_t floor;
  read_sm(from, offset, uint16_t, floor);

  uint32_t segments;
  read_sm(from, offset, uint32_t, segments);

  auto* indices = (uint32_t*) (from + offset);
  offset += segments * sizeof(uint32_t);

  auto* lengths = (uint32_t*) (from + offset);
  offset += segments * sizeof(uint32_t);

  return Stripped(total_l, actives, actives_l, indices, lengths, segments, floor);
}

uint32_t Packed::ZigZag(const int32_t x) {
  return (x << 1) ^ (x >> 31); // NOLINT(hicpp-signed-bitwise)
}

int32_t* Packed::DeltaEncode(const uint16_t* xs, const uint32_t len) {
  auto* diffed = static_cast<int32_t*>(std::malloc(len * sizeof(uint32_t)));
  std::adjacent_difference(xs, xs + len, diffed);
  return diffed;
}

uint8_t Packed::BestBits(const int32_t* xs, const uint32_t len) {
  std::array<size_t, 17> bit_counts = {}; // TODO replace with calloc?
  for (uint32_t i = 0; i < len; i++) {
    if (xs[i] == 0) { // log2(0) = -inf
      bit_counts[1]++;
    } else {
      bit_counts[std::log2(std::abs(xs[i])) + 1 + 1]++;
    }
  }

  uint8_t best_bits = 0;
  uint64_t best_size = std::numeric_limits<uint64_t>::max();
  for (uint32_t i = 1; i < 17; i++) {
    bit_counts[i] += bit_counts[i - 1];
    uint64_t size = i * bit_counts[i] + OUTLIER_BITS * (len - bit_counts[i]);
//    dbg("bestb: i=%" PRIu32 " size=%" PRIu64 "=%" PRIu32 "*%zu+%d*(%" PRIu32 "-%zu)", i, size, i, bit_counts[i], OUTLIER_BITS, len, bit_counts[i]);
    if (size < best_size) {
      best_bits = i;
      best_size = size;
    }
  }

  return best_bits;
}

uint64_t Packed::Write(char* to) const {
  uint64_t offset = 0;

  offset += signal.Write(to + offset);
  offset += outliers.Write(to + offset);

  return offset;
}

Packed::Packed(const Bitpacked &signal, const Bitpacked &outliers) : signal(signal), outliers(outliers) {}

Packed Packed::Pack(const uint16_t* xs, const uint32_t len) {
  dbg("pack:  len=%" PRIu32, len);
  int32_t* d_encoded = DeltaEncode(xs, len);

  uint8_t packed_bits = BestBits(d_encoded, len);

  const uint32_t outlier_marker = 0xffffffff >> (32u - packed_bits);
  // The marker is 00..0011...11 with bits_needed 1s
  // Would probably work with 11...11 too but I'm scared to try
  // TODO try ^
  // Everything below that top spot is fair game
  const uint32_t bound = outlier_marker - 1;
  dbg("pack:  bound=%" PRIu32, bound);

  // TODO replace with malloc & friends, perhaps
  std::vector<uint32_t> zigzagged;
  zigzagged.reserve(len); // TODO is this dumb?
  std::vector<uint32_t> outlier_list; // TODO initial size, too?
  for (uint32_t i = 0; i < len; i++) {
    uint32_t zigged = ZigZag(d_encoded[i]);
    dbg("pack:  zigged=%" PRIu32, zigged);
    if (zigged > bound) {
      outlier_list.push_back(zigged);
      zigzagged.push_back(outlier_marker);
    } else {
      zigzagged.push_back(zigged);
    }
  }

  dbg("pack:  len(zigzagged)=%zu len(outlier_list)=%zu", zigzagged.size(), outlier_list.size());

  Bitpacked signal = Bitpacked::Pack(zigzagged, packed_bits);
  Bitpacked outliers = Bitpacked::Pack(outlier_list, 17);

  return Packed(signal, outliers);
}

Packed Packed::Read(const char* from, uint64_t &offset) {
  Bitpacked signal = Bitpacked::Read(from, offset);
  Bitpacked outliers = Bitpacked::Read(from, offset);

  return Packed(signal, outliers);
}

// TODO I feel like the len and pointer are usually reversed?
uint16_t* Packed::Unpack(uint32_t* len) {
  uint32_t* u_outliers = outliers.Unpack();
  uint32_t* zigzagged = signal.Unpack();
  *len = signal.GetLen();
  auto* xs = static_cast<uint16_t*>(std::malloc(signal.GetLen() * sizeof(uint16_t))); // TODO don't love this

  uint32_t outlier_needle = 0;
  const uint32_t outlier_marker = 0xffffffff >> (32u - signal.GetBits());
  uint16_t prev = 0;
  for (uint32_t i = 0; i < signal.GetLen(); i++) {
    int32_t delta;
    if (zigzagged[i] == outlier_marker) {
      delta = UnZig(u_outliers[outlier_needle++]);
      dbg("upack: outlier=ye i=%" PRIu32 " prev=%" PRIu32 " delta=%" PRIi32, i, prev, delta);
    } else {
      delta = UnZig(zigzagged[i]);
      dbg("upack: outlier=no i=%" PRIu32 " prev=%" PRIu32 " delta=%" PRIi32, i, prev, delta);
    }
    xs[i] = prev + delta;
    dbg("upack: xs[%" PRIu32 "]=%" PRIu16, i, xs[i]);
    prev = xs[i];
  }

  return xs;
}

int32_t Packed::UnZig(const uint32_t x) {
  return (x >> 1u) ^ -(x & 1u);
}

uint64_t Bitpacked::Write(char* to) const {
  uint64_t offset = 0;

  dbg("write: len=%" PRIu32 " bits=%" PRIu8 " bytes=%" PRIu32, len, bits, bytes);

  write_sm(to, offset, uint32_t, len);
  write_sm(to, offset, uint8_t, bits);
  write_sm(to, offset, uint32_t, bytes);
  write_lg(to, offset, bytes, packed);

  return offset;
}

Bitpacked::Bitpacked(uint8_t* packed, uint32_t len, uint8_t bits, uint32_t bytes)
    : packed(packed), len(len), bits(bits), bytes(bytes) {}

Bitpacked Bitpacked::Pack(const std::vector<uint32_t> &xs, uint8_t bits) {
  const uint32_t pack_space = simdpack_compressedbytes(xs.size(), bits);
  auto* packed = new uint8_t[pack_space];
  __m128i* end_of_buf = simdpack_length(xs.data(), xs.size(), (__m128i*) packed, bits);
  uint32_t bytes = (end_of_buf - (__m128i*) packed) * sizeof(__m128i);
  return Bitpacked(packed, xs.size(), bits, bytes);
}

Bitpacked Bitpacked::Read(const char* from, uint64_t &offset) {
  uint32_t len;
  read_sm(from, offset, uint32_t, len);

  uint8_t bits;
  read_sm(from, offset, uint8_t, bits);

  uint32_t bytes;
  read_sm(from, offset, uint32_t, bytes);

  dbg("read:  len=%" PRIu32 " bits=%" PRIu8 " bytes=%" PRIu32, len, bits, bytes);

  auto* packed = (uint8_t*) (from + offset);
  offset += bytes;

  return Bitpacked(packed, len, bits, bytes);
}

uint32_t* Bitpacked::Unpack() {
  auto* unpacked = static_cast<uint32_t*>(malloc(len * sizeof(uint32_t)));
  simdunpack_length(reinterpret_cast<const __m128i*>(packed), len, unpacked, bits);
  return unpacked;
}

uint8_t Bitpacked::GetBits() const {
  return bits;
}

uint32_t Bitpacked::GetLen() const {
  return len;
}
