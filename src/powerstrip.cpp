// Just C dependencies
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>

// Used in debug mode
#include <cinttypes>
#include <cstdio>

#include "powerstrip.h"

#include "../extern/zstd/lib/common/huf.h" // Huffman coding
#include "../extern/simdcomp/include/simdcomp.h" // for bit-packing

#define WINDOW 3

#define OUTLIER_BITS 17

// A floor at this wattage and above won't be found
#define MAX_FLOOR 1000

#ifndef NDEBUG // cross-plat!
#define dbg(...) do {\
    std::printf("  ");\
    std::printf(__VA_ARGS__);\
    std::printf("\n");\
} while(0)
#else
#define dbg(...) do {} while(0);
#endif

struct stripped {
  uint32_t segments;
  uint32_t* indices;
  uint32_t* lengths;
  uint32_t actives_l;
  uint16_t floor;
  uint32_t total_l;
};

struct bitpacked {
  uint8_t* packed;
  uint32_t len;
  uint8_t bits;
  uint32_t bytes;
};

struct packed {
  struct bitpacked signal;
  struct bitpacked outliers;
};

bool is_active(const uint16_t x, const uint16_t floor) {
  return floor == UINT16_MAX || (x > floor + WINDOW || x < floor - WINDOW);
}

uint16_t find_floor(const uint16_t* xs, const uint32_t len) {
  auto* histogram = static_cast<uint32_t*>(std::calloc(MAX_FLOOR, sizeof(uint32_t)));
  for (uint32_t i = 0; i < len; i++) {
    if (xs[i] < MAX_FLOOR) {
      histogram[xs[i]]++;
    }
  }

  uint16_t floor = 0;
  uint32_t count = 0;
  for (uint32_t i = 0; i < MAX_FLOOR; i++) {
    if (histogram[i] > count) {
      floor = i;
      count = histogram[i];
    }
  }

  delete[] histogram;

  dbg("floor: floor=%" PRIu16 " count=%" PRIu32, floor, count);

  if (10 * count > len) {
    return floor;
  } else {
    return UINT16_MAX; // no floor with >10% of samples found
  }
}

struct stripped strip(const uint16_t* xs, const uint32_t len, uint32_t* actives) {
  struct stripped s{
      .indices = new uint32_t[len],
      .lengths = new uint32_t[len],
      .total_l = len,
  };

  uint16_t floor = find_floor(xs, len);

  bool was_active = false;
  uint32_t segments = 0;
  uint32_t segment_l = 0;
  uint32_t actives_l = 0;

  if (len >= 1) dbg("strip: xs[0]=%" PRIu16 " len=%" PRIu32, xs[0], len);

  for (size_t i = 0; i < len; i++) {
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

struct bitpacked bitpack(const uint32_t* xs, const uint32_t len, const uint8_t bits) {
  const uint32_t pack_space = simdpack_compressedbytes(len, bits);
  auto* packed = new uint8_t[pack_space];
  __m128i* end_of_buf = simdpack_length(xs, len, (__m128i*) packed, bits);
  uint32_t bytes = (end_of_buf - (__m128i*) packed) * sizeof(__m128i);
  struct bitpacked b{packed, len, bits, bytes};
  return b;
}

void delta_encode(uint32_t* xs, uint32_t len) {
  int32_t prev = 0;
  for (uint32_t i = 0; i < len; i++) {
    int32_t tmp_prev = xs[i];
    xs[i] -= prev;
    prev = tmp_prev;
  }
}

uint8_t best_bits(const int32_t* xs, const uint32_t len) {
  uint32_t bit_counts[17] = {0};

  for (uint32_t i = 0; i < len; i++) {
    if (xs[i] == 0) { // log2(0) = -inf
      bit_counts[1]++;
    } else {
      const uint8_t needed = std::log2(std::abs(xs[i])) + 1 + 1;
      bit_counts[needed]++;
    }
  }

  uint8_t best_bits = 0;
  uint64_t best_size = UINT64_MAX;
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
  const uint32_t bound = outlier_marker - 1;
  dbg("pack:  bound=%" PRIu32, bound);

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

  dbg("read:  len=%" PRIu32 " bits=%" PRIu8 " bytes=%" PRIu32, b->len, b->bits, b->bytes);

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

uint64_t read_stripped(const char* block, struct stripped* s) {
  uint64_t offset = 0;

  read_sm(block, offset, uint32_t, s->total_l);
  read_sm(block, offset, uint16_t, s->floor);
  read_sm(block, offset, uint32_t, s->segments);
  s->indices = (uint32_t*) (block + offset);
  offset += s->segments * sizeof(s->indices[0]);
  s->lengths = (uint32_t*) (block + offset);
  offset += s->segments * sizeof(s->lengths[0]);

  dbg("read:  total_l=%" PRIu32 " floor=%" PRIu16 " segments=%" PRIu32, s->total_l, s->floor, s->segments);

  return offset;
}

uint64_t uncompressed_size(const struct stripped* s, const struct packed* p) {
  return
      sizeof(p->signal.len) + sizeof(p->signal.bits) + sizeof(p->signal.bytes) + p->signal.bytes +
          sizeof(p->outliers.len) + sizeof(p->outliers.bits) + sizeof(p->outliers.bytes) + p->outliers.bytes +
          sizeof(s->total_l) + sizeof(s->floor) + sizeof(s->segments) +
          s->segments * sizeof(s->indices[0]) + s->segments * sizeof(s->lengths[0]);
}

uint64_t internal_compress(const uint16_t* block, const size_t len, char* out) {
  auto* actives = new uint32_t[len]; // some unfortunate special-casing because we need uint32_t, not uint16_t
  const struct stripped s = strip(block, len, actives);
  const struct packed p = pack(actives, s.actives_l);

  const uint64_t size = uncompressed_size(&s, &p);

  // If over INPUT_SIZE, we didn't do a good job. Just return plain uint16s
  if (size > INPUT_SIZE) {
    return 0;
  }

  uint64_t offset = 0;
  offset += write_bitpacked(&p.signal, out + offset);
  offset += write_bitpacked(&p.outliers, out + offset);
  offset += write_stripped(&s, out + offset);

  assert(offset == size);

  delete[] actives;
  delete[] s.lengths;
  delete[] s.indices;
  delete[] p.signal.packed;
  delete[] p.outliers.packed;

  return offset;
}

uint64_t write_tagged(const uint32_t tag, const char* xs, const size_t len, char* to) {
  uint64_t offset = 0;
  write_sm(to, offset, uint32_t, tag);
  std::memcpy(to + offset, xs, len);
  offset += len;
  return offset;
}

uint64_t compress_block(const uint16_t* block, const size_t len, char* out) {
  dbg("~");
  dbg("comp:  len=%zu", len);

  auto* uncompressed = new char[INPUT_SIZE];
  uint64_t uncompressed_size = internal_compress(block, len, uncompressed);

  if (uncompressed_size == 0) {
    dbg("compb: tag=0x00000000");
    uint64_t out_size = write_tagged(0x00000000, reinterpret_cast<const char*>(block), len * sizeof(uint16_t), out);
    delete[] uncompressed;
    return out_size;
  }

  size_t
      huf_size = HUF_compress(out + sizeof(uint32_t), OUTPUT_SIZE - sizeof(uint32_t), uncompressed, uncompressed_size);
  if (huf_size == 0 || HUF_isError(huf_size)) {
    dbg("HUF error: %s", HUF_isError(huf_size) ? HUF_getErrorName(huf_size) : "uncompressible");
    dbg("srcSize=%zu", (size_t) uncompressed_size);
    dbg("compb: tag=0xffffffff uncompressed_size=%" PRIu64, uncompressed_size);
    uint64_t out_size = write_tagged(0xffffffff, uncompressed, uncompressed_size, out);
    delete[] uncompressed;
    return out_size;
  }

  delete[] uncompressed;

  dbg("compb: tag=%" PRIu64, uncompressed_size);
  uint64_t offset = 0;
  write_sm(out, offset, uint32_t, (uint32_t) uncompressed_size);
  return sizeof(uint32_t) + huf_size;
}

uint64_t decompress_internal(const char* block, const size_t len, uint16_t* out) {
  struct packed p; // NOLINT(cppcoreguidelines-pro-type-member-init)
  struct stripped s; // NOLINT(cppcoreguidelines-pro-type-member-init)

  uint64_t offset = 0;
  offset += read_bitpacked(block + offset, &p.signal);
  offset += read_bitpacked(block + offset, &p.outliers);
  offset += read_stripped(block + offset, &s);
  assert(offset == len);

  auto* actives = new uint16_t[INPUT_SIZE / sizeof(uint16_t)];
  unpack(&p, actives);
  unstrip(&s, actives, out);
  delete[] actives;

  return s.total_l;
}

uint64_t decompress_block(const char* block, const size_t len, uint16_t* out) {
  uint64_t offset = 0;
  uint32_t tag;
  read_sm(block, offset, uint32_t, tag);

  if (tag == 0x00000000) { // plain uint16s
    dbg("decob: tag=0x00000000");
    std::memcpy(out, block + offset, len - offset);
    return (len - offset) / sizeof(uint16_t);
  } else if (tag == 0xffffffff) { // plain powerstripped
    dbg("decob: tag=0xffffffff");
    return decompress_internal(block + offset, len - offset, out);
  } else { // huffman coded and powerstripped
    dbg("decob: tag=%" PRIu32, tag);
    char* decomped = new char[INPUT_SIZE];
    size_t original_size = HUF_decompress(decomped, tag, block + offset, len - offset);
    assert(original_size == tag);
    uint64_t out_size = decompress_internal(decomped, original_size, out);
    delete[] decomped;
    return out_size;
  }
}
