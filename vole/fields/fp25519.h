#ifndef VOLE_FIELDS_FP25519_H__
#define VOLE_FIELDS_FP25519_H__

// Finite field F_p, p = 2^255 - 19 (the Curve25519 base field).
//
//   - Storage: four little-endian 64-bit limbs (`U256`), always reduced to
//              [0, p). One spare bit, but additions are still reduced at once
//              (lazy_adds = 1); the LPN encoders use the Acc below.
//   - Add/sub: a + b < 2p < 2^256, so a + b >= p iff a + b + 19 >= 2^255;
//              a borrow in a - b adds p back. Branch-free selects.
//   - Mul:     4x4 schoolbook (16 64x64->128 products) into eight limbs, then
//              2^256 = 38 (mod p): (L + 38 H), a small carry fold, and one
//              conditional subtract.
//   - Inverse: Fermat, a^(p-2).
//   - Sampling: 256 uniform bits masked to 255 bits are >= p with probability
//              19 / 2^255, so one conditional subtract is 2^-250-close to
//              uniform; a GGM leaf expands to two blocks (from_2blocks).
//
// A 128-bit block maps to its integer value (< 2^128): fine as a hash seed or
// a commitment coefficient with <= 128 bits of entropy, NOT a uniform element.

#include "emp-tool/emp-tool.h"
#include "vole/fields/utils.h"
#include "vole/fields/u256.h"
#include <cstdint>
#include <cstring>

using namespace emp;

class FP25519 {
 public:
  using u128 = unsigned __int128;

  static constexpr uint64_t PR_bit_len       = 255;
  static constexpr uint64_t slot_stride_bits = 255;   // COPE decomposes Delta < 2^255
  static constexpr U256     PR               = {{~0ULL - 18, ~0ULL, ~0ULL, ~0ULL >> 1}};
  // Only referenced by the packed-field (PR_num_pack > 1) branches, which are
  // dead for this type; kept so those branches still compile.
  static constexpr uint64_t PR_mask     = ~0ULL;
  static constexpr uint64_t PR_num_pack = 1;
  static constexpr uint64_t DS_byte_len = 32;
  static constexpr uint64_t TOP_MASK    = ~0ULL >> 1;   // low 63 bits of limb 3

  U256 val;

  FP25519() : val{{0, 0, 0, 0}} {}
  FP25519(const U256 &v) : val(v) { canon(val); }   // reduces any 256-bit value
  FP25519(uint64_t v) : val{{v, 0, 0, 0}} {}
  FP25519(int v) : val{{0, 0, 0, 0}} {
    if (v >= 0) val.w[0] = (uint64_t)v;
    else { val.w[0] = (uint64_t)(-(int64_t)v); val = neg(val); }
  }

  void operator=(const FP25519 &rhs) { val = rhs.val; }
  void operator=(uint64_t rhs) { val = U256{{rhs, 0, 0, 0}}; }

  void setZero() { val = U256{{0, 0, 0, 0}}; }

  // Whether this encoding is the canonical one (a peer may send any bytes).
  bool canonical() const { U256 t; return u256::sub(t, val, PR) != 0; }
  void assign_no_mod(const U256 &rhs) { val = rhs; canon(val); }
  void assign_no_mod(uint64_t rhs) { val = U256{{rhs, 0, 0, 0}}; }

  U256 value() const { return val; }

  // SPDZ layout helpers (value and MAC live in the pair type); identity here.
  void setHigh(const FP25519 &rhs) { val = rhs.val; }
  FP25519 getHigh() const { return *this; }
  void setLow(const FP25519 &rhs) { val = rhs.val; }
  FP25519 getLow() const { return *this; }

  bool operator==(const FP25519 &rhs) const { return val == rhs.val; }
  bool operator!=(const FP25519 &rhs) const { return val != rhs.val; }
  bool operator==(uint64_t rhs) const { return val == U256{{rhs, 0, 0, 0}}; }
  bool operator!=(uint64_t rhs) const { return !(*this == rhs); }

  FP25519 operator+(const FP25519 &rhs) const { FP25519 r; r.val = add_mod(val, rhs.val); return r; }
  FP25519 operator-(const FP25519 &rhs) const { FP25519 r; r.val = sub_mod(val, rhs.val); return r; }
  FP25519 operator-() const { FP25519 r; r.val = neg(val); return r; }
  FP25519 operator*(const FP25519 &rhs) const { FP25519 r; r.val = mult_mod(val, rhs.val); return r; }

  FP25519 operator+(uint64_t rhs) const { return *this + FP25519(rhs); }
  FP25519 operator-(uint64_t rhs) const { return *this - FP25519(rhs); }
  FP25519 operator*(uint64_t rhs) const { return *this * FP25519(rhs); }

  // Multiply by 2^n_pos, 0 <= n_pos < 256 (COPE's bit recomposition).
  FP25519 operator<<(const int n_pos) const {
    U256 t{{0, 0, 0, 0}};
    t.w[n_pos >> 6] = 1ULL << (n_pos & 63);
    FP25519 r; r.val = mult_mod(val, t); return r;
  }

  FP25519 &operator+=(const FP25519 &rhs) { *this = *this + rhs; return *this; }
  FP25519 &operator-=(const FP25519 &rhs) { *this = *this - rhs; return *this; }
  FP25519 &operator*=(const FP25519 &rhs) { *this = *this * rhs; return *this; }

  static constexpr unsigned lazy_adds = 1;
  void add_raw(const FP25519 &rhs) { val = add_mod(val, rhs.val); }
  void reduce() {}

  FP25519 negate() const { return -*this; }

  // a^(p-2); returns 0 for 0.
  FP25519 inv() const {
    const U256 e{{~0ULL - 20, ~0ULL, ~0ULL, ~0ULL >> 1}};  // p - 2
    U256 r{{1, 0, 0, 0}};
    for (int i = 254; i >= 0; --i) {
      r = mult_mod(r, r);
      if (e.bit(i)) r = mult_mod(r, val);
    }
    FP25519 o; o.val = r; return o;
  }

  template <typename PRNG>
  void rand(PRNG &prg) {
    prg.random_data_unaligned(&val, sizeof(U256));
    from_255(val);
  }
  // n uniform elements from one PRG call (COPE's per-bit streams).
  template <typename PRNG>
  static void prg_elements(PRNG &prg, FP25519 *out, std::size_t n) {
    prg.random_data_unaligned(out, n * sizeof(U256));
    for (std::size_t i = 0; i < n; ++i) from_255(out[i].val);
  }

  template <typename IO> void send(IO *netio) { netio->send_data(&val, sizeof(U256)); }
  template <typename IO> void recv(IO *netio) { netio->recv_data(&val, sizeof(U256)); canon(val); }

  block hash() { return Hash::hash_for_block(&val, sizeof(U256)); }

  void from_block(block rhs) {
    std::memcpy(&val.w[0], &rhs, sizeof(block));
    val.w[2] = val.w[3] = 0;
  }
  // 255 uniform bits from two blocks: uniform up to 2^-250.
  void from_2blocks(block b0, block b1) {
    std::memcpy(&val.w[0], &b0, sizeof(block));
    std::memcpy(&val.w[2], &b1, sizeof(block));
    from_255(val);
  }

  static U256 copy_compose(const U256 &basis) { return basis; }
  static uint64_t copy_compose(uint64_t basis) { return basis; }

  // ---- limb arithmetic -----------------------------------------------------

  // x < 2p -> x mod p: x >= p iff x + 19 >= 2^255 (x + 19 < 2^256, and for
  // x >= p clearing bit 255 of x + 19 gives x - p).
  static void canon255(U256 &x) {
    U256 t = x;
    u128 c = (u128)t.w[0] + 19;
    t.w[0] = (uint64_t)c;
    for (int i = 1; i < 4; ++i) { c = (u128)t.w[i] + (uint64_t)(c >> 64); t.w[i] = (uint64_t)c; }
    uint64_t ge = 0 - (t.w[3] >> 63);
    t.w[3] &= TOP_MASK;
    u256::select(x, t, ge);
  }
  // Uniform 256 bits -> uniform element: drop bit 255, then canon255.
  static void from_255(U256 &x) {
    x.w[3] &= TOP_MASK;
    canon255(x);
  }
  // Any x < 2^256 -> x mod p: fold bit 255 (2^255 = 19), then canon255.
  static void canon(U256 &x) {
    uint64_t hi = x.w[3] >> 63;
    x.w[3] &= TOP_MASK;
    u128 c = (u128)x.w[0] + 19 * hi;
    x.w[0] = (uint64_t)c;
    for (int i = 1; i < 4; ++i) { c = (u128)x.w[i] + (uint64_t)(c >> 64); x.w[i] = (uint64_t)c; }
    canon255(x);  // x < 2^255 + 19 here; canon255 also handles bit 255 set
  }
  // x + k * 2^256 (k < 2^57) -> x mod p, with 2^256 = 38.
  static void fold(U256 &x, uint64_t k) {
    u128 c = (u128)x.w[0] + (u128)38 * k;
    x.w[0] = (uint64_t)c;
    for (int i = 1; i < 4; ++i) { c = (u128)x.w[i] + (uint64_t)(c >> 64); x.w[i] = (uint64_t)c; }
    uint64_t c2 = (uint64_t)(c >> 64);  // 0 or 1; if 1, x < 38 * (k + 1) now
    c = (u128)x.w[0] + 38 * c2;
    x.w[0] = (uint64_t)c;
    for (int i = 1; i < 4; ++i) { c = (u128)x.w[i] + (uint64_t)(c >> 64); x.w[i] = (uint64_t)c; }
    canon(x);
  }

  static U256 add_mod(const U256 &a, const U256 &b) {
    U256 s;
    u256::add(s, a, b);  // < 2p < 2^256: no carry out
    canon255(s);         // s < 2^256 with s < 2p: s >= p iff s + 19 >= 2^255
    return s;
  }

  static U256 sub_mod(const U256 &a, const U256 &b) {
    U256 d, t;
    uint64_t borrow = u256::sub(d, a, b);
    u256::add(t, d, PR);
    u256::select(d, t, 0 - borrow);
    return d;
  }

  static U256 neg(const U256 &a) { return sub_mod(U256{{0, 0, 0, 0}}, a); }

  static U256 mult_mod(const U256 &a, const U256 &b) {
    uint64_t c[8];
    u256::mul_wide(c, a, b);
    // r = L + 38 * H, L = c[0..3], H = c[4..7]
    U256 r;
    u128 t = 0;
    for (int i = 0; i < 4; ++i) {
      t = (u128)c[i] + (u128)38 * c[4 + i] + (uint64_t)(t >> 64);
      r.w[i] = (uint64_t)t;
    }
    fold(r, (uint64_t)(t >> 64));  // top < 39
    return r;
  }

  static std::size_t size() { return sizeof(U256); }

  // Lazy accumulator for sums of many elements (the LPN encoders): one
  // 128-bit running sum per limb, carries resolved and folded once in get().
  struct Acc {
    u128 s[4];
    explicit Acc(const FP25519 &x) : s{x.val.w[0], x.val.w[1], x.val.w[2], x.val.w[3]} {}
    void add(const FP25519 &y) {
      for (int i = 0; i < 4; ++i) s[i] += y.val.w[i];
    }
    FP25519 get() const {
      U256 r;
      u128 t = s[0];
      r.w[0] = (uint64_t)t;
      for (int i = 1; i < 4; ++i) { t = s[i] + (t >> 64); r.w[i] = (uint64_t)t; }
      FP25519 o;
      o.val = r;
      fold(o.val, (uint64_t)(t >> 64));
      return o;
    }
  };
};
static_assert(sizeof(FP25519) == 32, "FP25519 must be 32 bytes (sent/hashed as raw bytes)");

// GGM leaves expand to two blocks (spfss_*.h).
template <> struct wide_leaf<FP25519> { static constexpr bool value = true; };
// Commitment-matrix coefficients come from from_block on a PRP block whose top
// 32 bits also pick the row (com_matrix.h): >= 96 uniform bits per row.
template <> struct com_coef_bits<FP25519> { static constexpr double value = 96.0; };
// MPFSS check: coefficient-array form (long multiplication, as FP192).
template <> struct poly_hash_horner<FP25519> { static constexpr bool value = false; };

#endif  // VOLE_FIELDS_FP25519_H__
