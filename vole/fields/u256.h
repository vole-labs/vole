#ifndef VOLE_FIELDS_U256_H__
#define VOLE_FIELDS_U256_H__

// Raw 256-bit integer, four little-endian 64-bit limbs: the storage of the
// 256-bit prime fields (fp25519.h, fpp256.h). Plain aggregate, 32 bytes, so
// arrays of it can be filled by a PRG and sent or hashed as bytes.

#include <cstddef>
#include <cstdint>

struct U256 {
  uint64_t w[4];
  bool operator==(const U256 &o) const {
    return w[0] == o.w[0] && w[1] == o.w[1] && w[2] == o.w[2] && w[3] == o.w[3];
  }
  bool operator!=(const U256 &o) const { return !(*this == o); }
  bool bit(std::size_t i) const { return (w[i >> 6] >> (i & 63)) & 1; }
};
static_assert(sizeof(U256) == 32, "U256 must be 32 bytes");

namespace u256 {
using u128 = unsigned __int128;

// x = mask ? t : x, without a branch.
inline void select(U256 &x, const U256 &t, uint64_t mask) {
  for (int i = 0; i < 4; ++i) x.w[i] = (t.w[i] & mask) | (x.w[i] & ~mask);
}

// r = a + b mod 2^256; returns the carry out.
inline uint64_t add(U256 &r, const U256 &a, const U256 &b) {
  u128 t = 0;
  for (int i = 0; i < 4; ++i) {
    t = (u128)a.w[i] + b.w[i] + (uint64_t)(t >> 64);
    r.w[i] = (uint64_t)t;
  }
  return (uint64_t)(t >> 64);
}

// r = a - b mod 2^256; returns the borrow out (0 or 1).
inline uint64_t sub(U256 &r, const U256 &a, const U256 &b) {
  uint64_t borrow = 0;
  for (int i = 0; i < 4; ++i) {
    u128 t = (u128)a.w[i] - b.w[i] - borrow;
    r.w[i] = (uint64_t)t;
    borrow = (uint64_t)(t >> 64) & 1;
  }
  return borrow;
}

// 4x4 schoolbook product into eight limbs.
inline void mul_wide(uint64_t c[8], const U256 &a, const U256 &b) {
  for (int i = 0; i < 8; ++i) c[i] = 0;
  for (int i = 0; i < 4; ++i) {
    uint64_t carry = 0;
    for (int j = 0; j < 4; ++j) {
      u128 t = (u128)a.w[i] * b.w[j] + c[i + j] + carry;
      c[i + j] = (uint64_t)t;
      carry = (uint64_t)(t >> 64);
    }
    c[i + 4] = carry;
  }
}
}  // namespace u256

#endif  // VOLE_FIELDS_U256_H__
