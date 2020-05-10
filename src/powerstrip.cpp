#include <cstring> // for memcpy
#include <cinttypes> // for printf macros?
#include <numeric>
#include <array>
#include <limits>
#include <common.h>

#include "powerstrip.h"

#include "../extern/zstd/lib/common/huf.h" // Huffman coding
#include "../extern/simdcomp/include/simdcomp.h" // for bit-packing

#define WINDOW 0

#define OUTLIER_BITS 17

bool is_active(const uint16_t x, const uint16_t floor) {
  return x > floor + WINDOW || x < floor - WINDOW;
}

struct stripped strip(const uint16_t* xs, const uint32_t len, uint32_t* actives) {
  struct stripped s{
      .indices = new uint32_t[len],
      .lengths = new uint32_t[len],
      .total_l = len,
  };

  uint16_t floor = 0; // TODO do this

  bool was_active = false;
  uint32_t segments = 0;
  uint32_t segment_l = 0;
  uint32_t actives_l = 0;

  if (len >= 1) dbg("strip: xs[0]=%" PRIu16 " len=%" PRIu32, xs[0], len);

  for (size_t i = 0; i < len; i++) {
    // TODO consider inverting the nesting so it's more like an FSM
    if (is_active(xs[i], floor)) {
      if (!was_active) {
        segment_l = 0;
        s.indices[segments++] = i;
        was_active = true;
      }
      actives[actives_l++] = xs[i];
      segment_l++;
    } else if (was_active) {
      s.lengths[segments - 1] = segment_l;
      was_active = false;
    }
  }

  if (was_active) {
    s.lengths[segments - 1] = segment_l;
  }

  dbg("strip: actives_l=%" PRIu32 " segments=%" PRIu32, actives_l, segments);

  s.segments = segments;
  s.actives_l = actives_l;
  s.floor = floor;

  return s;
}

void unstrip(const struct stripped* s, const uint16_t* actives, uint16_t* to) {
  for (uint32_t i = 0; i < s->total_l; i++) {
    to[i] = s->floor;
  }

  uint32_t needle = 0;
  for (uint32_t i = 0; i < s->segments; i++) {
    for (uint32_t j = 0; j < s->lengths[i]; j++) {
      dbg("ustrp: to[indices[%" PRIu32 "]+%" PRIu32 "]=actives[%" PRIu32 "]=%" PRIu16, i, j, needle, actives[needle]);
      to[s->indices[i] + j] = actives[needle++];
    }
  }
}

struct bitpacked bitpack(const uint32_t* xs, const uint32_t len, uint8_t bits) {
  const uint32_t pack_space = simdpack_compressedbytes(len, bits);
  auto* packed = new uint8_t[pack_space];
  __m128i* end_of_buf = simdpack_length(xs, len, (__m128i*) packed, bits);
  uint32_t bytes = (end_of_buf - (__m128i*) packed) * sizeof(__m128i);
  struct bitpacked b{packed, len, bits, bytes};
  return b;
}

void delta_encode(uint32_t* xs, uint32_t len) {
  std::adjacent_difference(xs, xs + len, xs);
}

uint8_t best_bits(const int32_t* xs, const uint32_t len) {
  uint32_t bit_counts[17] = {0};

  for (uint32_t i = 0; i < len; i++) {
    if (xs[i] == 0) { // log2(0) = -inf
      bit_counts[1]++;
    } else {
      // TODO log2 is slow--we can do a bitscan like simdcomp if we do zigzagging earlier
      const uint8_t needed = std::log2(std::abs(xs[i])) + 1 + 1;
      bit_counts[needed]++;
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

uint32_t zig_zag(const int32_t x) {
  return (x << 1) ^ (x >> 31); // NOLINT(hicpp-signed-bitwise)
}

int32_t un_zig(const uint32_t x) {
  return (x >> 1u) ^ -(x & 1u);
}

void unpack(const struct packed* p, uint16_t* actives) {
  auto* outliers = new uint32_t[p->outliers.len];
  simdunpack_length(reinterpret_cast<const __m128i*>(p->outliers.packed), p->outliers.len, outliers, p->outliers.bits);

  // TODO would be nice to just unpack into actives but uint32_t != uint16_t :l
  auto* zigged = new uint32_t[p->signal.len];
  simdunpack_length(reinterpret_cast<const __m128i*>(p->signal.packed), p->signal.len, zigged, p->signal.bits);

  uint32_t outlier_needle = 0;
  const uint32_t outlier_marker = 0xffffffff >> (32u - p->signal.bits);
  uint16_t prev = 0;
  for (uint32_t i = 0; i < p->signal.len; i++) {
    int32_t delta;
    if (zigged[i] == outlier_marker) {
      delta = un_zig(outliers[outlier_needle++]);
      dbg("upack: outlier=ye i=%" PRIu32 " prev=%" PRIu32 " delta=%" PRIi32, i, prev, delta);
    } else {
      delta = un_zig(zigged[i]);
      dbg("upack: outlier=no i=%" PRIu32 " prev=%" PRIu32 " delta=%" PRIi32, i, prev, delta);
    }
    actives[i] = prev + delta;
    dbg("upack: xs[%" PRIu32 "]=%" PRIu16, i, actives[i]);
    prev = actives[i];
  }

  delete[] zigged;
  delete[] outliers;
}

struct packed pack(uint32_t* xs, const uint32_t len) {
  dbg("pack:  len=%" PRIu32, len);
  delta_encode(xs, len);
  auto* d_encoded = (int32_t*) xs;

  uint8_t packed_bits = best_bits(d_encoded, len);

  const uint32_t outlier_marker = 0xffffffff >> (32u - packed_bits);
  // The marker is 00..0011...11 with bits_needed 1s
  // Would probably work with 11...11 too but I'm scared to try
  // TODO try ^ but this is not going to give you more than 1 cycle :l
  // Everything below that top spot is fair game
  const uint32_t bound = outlier_marker - 1;
  dbg("pack:  bound=%" PRIu32, bound);

  // TODO can use xs instead of zigzagged if we are grasping for memory
  auto* zigzagged = new uint32_t[len];
  auto* outliers = new uint32_t[len];
  uint32_t outliers_l = 0;
  for (uint32_t i = 0; i < len; i++) {
    uint32_t zigged = zig_zag(d_encoded[i]);
    if (zigged > bound) {
      outliers[outliers_l++] = zigged;
      zigzagged[i] = outlier_marker;
    } else {
      zigzagged[i] = zigged;
    }
  }

  dbg("pack:  len(zigzagged)=%" PRIu32 " len(outlier_list)=%" PRIu32, len, outliers_l);

  struct packed p{
      .signal = bitpack(zigzagged, len, packed_bits),
      .outliers = bitpack(outliers, outliers_l, 17),
  };

  delete[] zigzagged;
  delete[] outliers;

  return p;
}

#define write_sm(buf, offset, type, x) do {\
  (reinterpret_cast<type*>(buf+offset))[0] = x;\
  offset += sizeof(type);\
} while(0)

#define write_lg(buf, offset, size, x) do {\
  memcpy(buf + offset, x, size);\
  offset += size;\
} while (0)

#define read_sm(buf, offset, type, x) do {\
  x = (reinterpret_cast<const type*>(buf + offset))[0];\
  offset += sizeof(type);\
} while(0)

// TODO maybe pass by value instead?
uint64_t write_bitpacked(const struct bitpacked* b, char* to) {
  uint64_t offset = 0;

  dbg("write: len=%" PRIu32 " bits=%" PRIu8 " bytes=%" PRIu32, b->len, b->bits, b->bytes);

  write_sm(to, offset, uint32_t, b->len);
  write_sm(to, offset, uint8_t, b->bits);
  write_sm(to, offset, uint32_t, b->bytes);
  write_lg(to, offset, b->bytes, b->packed);

  return offset;
}

uint64_t read_bitpacked(const char* block, struct bitpacked* b) {
  uint64_t offset = 0;

  read_sm(block, offset, uint32_t, b->len);
  read_sm(block, offset, uint8_t, b->bits);
  read_sm(block, offset, uint32_t, b->bytes);
  b->packed = (uint8_t*) (block + offset);
  offset += b->bytes;

  return offset;
}

uint64_t write_stripped(const struct stripped* s, char* to) {
  uint64_t offset = 0;

  dbg("write: total_l=%" PRIu32 " floor=%" PRIu16 " segments=%" PRIu32, s->total_l, s->floor, s->segments);
  write_sm(to, offset, uint32_t, s->total_l);
  write_sm(to, offset, uint16_t, s->floor);
  write_sm(to, offset, uint32_t, s->segments);
  write_lg(to, offset, s->segments * sizeof(s->indices[0]), s->indices);
  write_lg(to, offset, s->segments * sizeof(s->lengths[0]), s->lengths);

  return offset;
}

uint16_t read_stripped(const char* block, struct stripped* s) {
  uint64_t offset = 0;

  read_sm(block, offset, uint32_t, s->total_l);
  read_sm(block, offset, uint16_t, s->floor);
  read_sm(block, offset, uint32_t, s->segments);
  s->indices = (uint32_t*) (block + offset);
  offset += s->segments * sizeof(s->indices[0]);
  s->lengths = (uint32_t*) (block + offset);
  offset += s->segments * sizeof(s->lengths[0]);

  return offset;
}

// TODO Should probably be uint32_t instead of size_t
uint64_t compress_block(const uint16_t* block, const size_t len, char* out) {
  dbg("~");
  dbg("comp:  len=%zu", len);

  auto* actives = new uint32_t[len]; // some unfortunate special-casing because we need uint32_t, not uint16_t
  struct stripped s = strip(block, len, actives);
  struct packed p = pack(actives, s.actives_l);

  uint64_t offset = 0;
  offset += write_bitpacked(&p.signal, out + offset);
  offset += write_bitpacked(&p.outliers, out + offset);
  offset += write_stripped(&s, out + offset);

  delete[] actives;
  delete[] s.lengths;
  delete[] s.indices;
  delete[] p.signal.packed;
  delete[] p.outliers.packed;

  return offset;
}

uint64_t decompress_block(const char* block, const size_t len, uint16_t* out) {
  // TODO initialize these or no? probs no
  struct packed p;
  struct stripped s;

  uint64_t offset = 0;
  offset += read_bitpacked(block + offset, &p.signal);
  offset += read_bitpacked(block + offset, &p.outliers);
  offset += read_stripped(block + offset, &s);
  assert(offset == len); // We read everything yay

  // TODO uneven len has possibility to break this
  // TODO or try BLOCK_SIZE? Can we make all mallocs of size BLOCK_SIZE?
  auto* actives = new uint16_t[BLOCK_SIZE / sizeof(uint16_t)];
  unpack(&p, actives);

  unstrip(&s, actives, out);

  delete[] actives;

  return s.total_l;
}
