#ifndef VOLE_FIELDS_FP128_H__
#define VOLE_FIELDS_FP128_H__

// Finite field F_p, p = 2^128 - 159 (the largest prime below 2^128).
//
//   - Storage: two little-endian 64-bit limbs (`U128`), always reduced to
//              [0, p). No spare bit: a sum of two elements can pass 2^128, so
//              additions are reduced at once (lazy_adds = 1); the LPN encoders
//              use the Acc below. Limb (class) storage also routes COPE and the
//              base sVOLE through their limb paths, as for FP192 / FP25519.
//   - Reduction: 2^128 = 159 (mod p). Any x < 2p satisfies x >= p iff
//              x + 159 carries out of 128 bits, and then (x + 159) mod 2^128
//              = x - p. Branch-free selects throughout.
//   - Mul:     2x2 schoolbook (four 64x64->128 products) into four limbs
//              L + 2^128 H, then L + 159 H (< 2^137), one more fold of the top
//              word by 159, and one conditional subtract.
//   - Inverse: Fermat, a^(p-2).
//   - Sampling: 128 uniform bits, one conditional subtract. A value >= p
//              occurs with probability 159 / 2^128 ~ 2^-120.7, which is the
//              statistical distance from uniform per element. This holds for
//              rand(), the COPE PRG outputs, and the GGM leaves (from_block on
//              one 128-bit block, so no wide leaves are needed).
//
// A 128-bit block maps to its integer value mod p (from_block), which is the
// near-uniform element above when the block is uniform.

#include "emp-tool/emp-tool.h"
#include "vole/fields/utils.h"
#include <cstdint>
#include <cstring>

using namespace emp;

// Raw 128-bit integer, two little-endian 64-bit limbs. Plain aggregate, 16
// bytes, so arrays of it can be filled by a PRG and sent or hashed as bytes.
struct U128 {
  uint64_t w[2];
  bool operator==(const U128 &o) const { return w[0] == o.w[0] && w[1] == o.w[1]; }
  bool operator!=(const U128 &o) const { return !(*this == o); }
  bool bit(std::size_t i) const { return (w[i >> 6] >> (i & 63)) & 1; }
};
static_assert(sizeof(U128) == 16, "U128 must be 16 bytes");

class FP128 {
 public:
  using u128 = unsigned __int128;

  static constexpr uint64_t PR_bit_len       = 128;
  static constexpr uint64_t slot_stride_bits = 128;   // COPE decomposes Delta < 2^128
  static constexpr uint64_t C                = 159;   // 2^128 - p
  static constexpr U128     PR               = {{0ULL - C, ~0ULL}};
  // Only referenced by the packed-field (PR_num_pack > 1) branches, which are
  // dead for this type; kept so those branches still compile.
  static constexpr uint64_t PR_mask     = ~0ULL;
  static constexpr uint64_t PR_num_pack = 1;
  static constexpr uint64_t DS_byte_len = 16;

  U128 val;

  FP128() : val{{0, 0}} {}
  FP128(const U128 &v) : val(v) { canon(val); }   // reduces any 128-bit value
  FP128(uint64_t v) : val{{v, 0}} {}
  FP128(int v) : val{{0, 0}} {
    if (v >= 0) val.w[0] = (uint64_t)v;
    else { val.w[0] = (uint64_t)(-(int64_t)v); val = neg(val); }
  }

  void operator=(const FP128 &rhs) { val = rhs.val; }
  void operator=(uint64_t rhs) { val = U128{{rhs, 0}}; }

  void setZero() { val = U128{{0, 0}}; }

  // Whether this encoding is the canonical one (a peer may send any bytes).
  bool canonical() const { return !(val.w[1] == ~0ULL && val.w[0] >= PR.w[0]); }
  void assign_no_mod(const U128 &rhs) { val = rhs; canon(val); }
  void assign_no_mod(uint64_t rhs) { val = U128{{rhs, 0}}; }

  U128 value() const { return val; }

  // SPDZ layout helpers (value and MAC live in the pair type); identity here.
  void setHigh(const FP128 &rhs) { val = rhs.val; }
  FP128 getHigh() const { return *this; }
  void setLow(const FP128 &rhs) { val = rhs.val; }
  FP128 getLow() const { return *this; }

  bool operator==(const FP128 &rhs) const { return val == rhs.val; }
  bool operator!=(const FP128 &rhs) const { return val != rhs.val; }
  bool operator==(uint64_t rhs) const { return val == U128{{rhs, 0}}; }
  bool operator!=(uint64_t rhs) const { return !(*this == rhs); }

  FP128 operator+(const FP128 &rhs) const { FP128 r; r.val = add_mod(val, rhs.val); return r; }
  FP128 operator-(const FP128 &rhs) const { FP128 r; r.val = sub_mod(val, rhs.val); return r; }
  FP128 operator-() const { FP128 r; r.val = neg(val); return r; }
  FP128 operator*(const FP128 &rhs) const { FP128 r; r.val = mult_mod(val, rhs.val); return r; }

  FP128 operator+(uint64_t rhs) const { return *this + FP128(rhs); }
  FP128 operator-(uint64_t rhs) const { return *this - FP128(rhs); }
  FP128 operator*(uint64_t rhs) const { return *this * FP128(rhs); }

  // Multiply by 2^n_pos, 0 <= n_pos < 128 (COPE's bit recomposition).
  FP128 operator<<(const int n_pos) const {
    U128 t{{0, 0}};
    t.w[n_pos >> 6] = 1ULL << (n_pos & 63);
    FP128 r; r.val = mult_mod(val, t); return r;
  }

  FP128 &operator+=(const FP128 &rhs) { *this = *this + rhs; return *this; }
  FP128 &operator-=(const FP128 &rhs) { *this = *this - rhs; return *this; }
  FP128 &operator*=(const FP128 &rhs) { *this = *this * rhs; return *this; }

  static constexpr unsigned lazy_adds = 1;
  void add_raw(const FP128 &rhs) { val = add_mod(val, rhs.val); }
  void reduce() {}

  FP128 negate() const { return -*this; }

  // a^(p-2); returns 0 for 0.
  FP128 inv() const {
    const U128 e{{0ULL - C - 2, ~0ULL}};  // p - 2
    U128 r{{1, 0}};
    for (int i = 127; i >= 0; --i) {
      r = mult_mod(r, r);
      if (e.bit(i)) r = mult_mod(r, val);
    }
    FP128 o; o.val = r; return o;
  }

  template <typename PRNG>
  void rand(PRNG &prg) {
    prg.random_data_unaligned(&val, sizeof(U128));
    canon(val);
  }

  template <typename IO> void send(IO *netio) { netio->send_data(&val, sizeof(U128)); }
  template <typename IO> void recv(IO *netio) { netio->recv_data(&val, sizeof(U128)); canon(val); }

  block hash() { return Hash::hash_for_block(&val, sizeof(U128)); }

  void from_block(block rhs) {
    std::memcpy(&val, &rhs, sizeof(block));
    canon(val);
  }

  static U128 copy_compose(const U128 &basis) { return basis; }
  static uint64_t copy_compose(uint64_t basis) { return basis; }

  // ---- limb arithmetic -----------------------------------------------------

  static u128 to(const U128 &x) { return ((u128)x.w[1] << 64) | x.w[0]; }
  static U128 of(u128 x) { return U128{{(uint64_t)x, (uint64_t)(x >> 64)}}; }

  // x = mask ? t : x, without a branch.
  static u128 select(u128 x, u128 t, uint64_t mask) {
    u128 m = ((u128)mask << 64) | mask;
    return (t & m) | (x & ~m);
  }

  // x < 2^128 (any value, hence < 2p) -> x mod p.
  static u128 canon128(u128 x) {
    u128 t = x + C;
    uint64_t ge = 0 - (uint64_t)(t < x);   // carry out: x >= p
    return select(x, t, ge);
  }
  static void canon(U128 &x) { x = of(canon128(to(x))); }

  // lo + k * 2^128 (k < 2^56) -> mod p, with 2^128 = 159.
  static u128 fold(u128 lo, uint64_t k) {
    u128 s = lo + (u128)(C * k);                // C * k < 2^64
    uint64_t c = (uint64_t)(s < lo);            // one carry out at most
    u128 s2 = s + (u128)(C * c);                // if c, s < C * k: no carry
    return canon128(s2);
  }

  static U128 add_mod(const U128 &a, const U128 &b) {
    u128 x = to(a), s = x + to(b);
    uint64_t c = (uint64_t)(s < x);             // a + b >= 2^128
    u128 t = s + C;
    uint64_t c2 = (uint64_t)(t < s);            // a + b in [p, 2^128)
    // c: true sum is s + 2^128 = s + 159 (mod p) and s + 159 < 2^128;
    // c2: s >= p and t = s - p.
    return of(select(s, t, 0 - (c | c2)));
  }

  static U128 sub_mod(const U128 &a, const U128 &b) {
    u128 x = to(a), y = to(b), d = x - y;
    uint64_t borrow = 0 - (uint64_t)(x < y);
    // borrow: d = a - b + 2^128, and a - b + p = d - 159 (>= 1, no wrap)
    return of(select(d, d - C, borrow));
  }

  static U128 neg(const U128 &a) { return sub_mod(U128{{0, 0}}, a); }

  static U128 mult_mod(const U128 &a, const U128 &b) {
    u128 p00 = (u128)a.w[0] * b.w[0];
    u128 p01 = (u128)a.w[0] * b.w[1];
    u128 p10 = (u128)a.w[1] * b.w[0];
    u128 p11 = (u128)a.w[1] * b.w[1];
    // c = c3 c2 c1 c0 = a * b
    uint64_t c0 = (uint64_t)p00;
    u128 t = (p00 >> 64) + (uint64_t)p01 + (uint64_t)p10;
    uint64_t c1 = (uint64_t)t;
    t = (t >> 64) + (p01 >> 64) + (p10 >> 64) + (uint64_t)p11;
    uint64_t c2 = (uint64_t)t;
    uint64_t c3 = (uint64_t)((t >> 64) + (p11 >> 64));
    // L + 159 H, L = (c1, c0), H = (c3, c2): < 2^137, top word < 2^9
    u128 m0 = (u128)C * c2, m1 = (u128)C * c3;
    t = (u128)c0 + (uint64_t)m0;
    uint64_t r0 = (uint64_t)t;
    t = (t >> 64) + c1 + (m0 >> 64) + (uint64_t)m1;
    uint64_t r1 = (uint64_t)t;
    uint64_t top = (uint64_t)((t >> 64) + (m1 >> 64));
    return of(fold(((u128)r1 << 64) | r0, top));
  }

  static std::size_t size() { return sizeof(U128); }

  // Lazy accumulator for sums of many elements (the LPN encoders): one
  // 128-bit running sum per limb, carries resolved and folded once in get().
  // Valid for fewer than 2^55 terms.
  struct Acc {
    u128 s[2];
    explicit Acc(const FP128 &x) : s{x.val.w[0], x.val.w[1]} {}
    void add(const FP128 &y) { s[0] += y.val.w[0]; s[1] += y.val.w[1]; }
    FP128 get() const {
      u128 hi = s[1] + (s[0] >> 64);
      u128 lo = ((u128)(uint64_t)hi << 64) | (uint64_t)s[0];
      FP128 o;
      o.val = of(fold(lo, (uint64_t)(hi >> 64)));
      return o;
    }
  };
};
static_assert(sizeof(FP128) == 16, "FP128 must be 16 bytes (sent/hashed as raw bytes)");

// Commitment-matrix coefficients come from from_block on a PRP block whose top
// 32 bits also pick the row (com_matrix.h): >= 96 uniform bits per row.
template <> struct com_coef_bits<FP128> { static constexpr double value = 96.0; };

#endif  // VOLE_FIELDS_FP128_H__
