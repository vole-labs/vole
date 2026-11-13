#ifndef EMP_UTILS_FIELD_FP192_H__
#define EMP_UTILS_FIELD_FP192_H__

// Finite field F_p, p = 2^192 - 2^64 - 1 (the NIST P-192 prime).
//
// Arithmetic backbone:
//   - Storage: three little-endian 64-bit limbs (`U192`), always reduced to
//              [0, p). There is no spare bit, so additions are reduced
//              immediately (lazy_adds = 1).
//   - Add/sub: carry/borrow chains; an overflow past 2^192 folds back as
//              2^192 = 2^64 + 1 (mod p).
//   - Mul:     3x3 schoolbook (9 64x64->128 products) into six limbs, then
//              the Solinas reduction of FIPS 186-4 D.2.1:
//                c = (c5..c0) -> (c2,c1,c0) + (0,c3,c3) + (c4,c4,0) + (c5,c5,c5)
//              followed by a small carry fold and one conditional subtract.
//   - Inverse: Fermat, a^(p-2).
//   - Sampling: 192 uniform bits are >= p with probability 2^-128, so one
//              conditional subtract gives a distribution 2^-128-close to
//              uniform.
//
// A 128-bit block maps to a field element as its integer value (< 2^128), so it
// is NOT a uniform element; the GGM leaves are therefore expanded to two blocks
// (`from_2blocks`, gated in spfss_* by `wide_leaf<FP>`).

#include "emp-tool/emp-tool.h"
#include "vole/fields/utils.h"
#include <cstdint>
#include <cstring>

using namespace emp;

// Raw 192-bit integer (little-endian limbs). Plain aggregate, 24 bytes, so
// arrays of it can be filled by a PRG and sent as bytes.
struct U192 {
  uint64_t w[3];
  bool operator==(const U192 &o) const {
    return w[0] == o.w[0] && w[1] == o.w[1] && w[2] == o.w[2];
  }
  bool operator!=(const U192 &o) const { return !(*this == o); }
  bool bit(std::size_t i) const { return (w[i >> 6] >> (i & 63)) & 1; }
};
static_assert(sizeof(U192) == 24, "U192 must be 24 bytes");

class FP192 {
 public:
  using u128 = unsigned __int128;

  static constexpr uint64_t PR_bit_len       = 192;
  static constexpr uint64_t slot_stride_bits = 192;   // = PR_bit_len
  static constexpr U192     PR               = {{~0ULL, ~0ULL - 1, ~0ULL}};
  // Only referenced by the packed-field (PR_num_pack > 1) branches, which are
  // dead for this type; kept so those branches still compile.
  static constexpr uint64_t PR_mask     = ~0ULL;
  static constexpr uint64_t PR_num_pack = 1;
  static constexpr uint64_t DS_byte_len = 24;

  U192 val;

  FP192() : val{{0, 0, 0}} {}
  FP192(const U192 &v) : val(v) { canon(val); }        // reduces
  FP192(uint64_t v) : val{{v, 0, 0}} {}
  FP192(int v) : val{{0, 0, 0}} {
    if (v >= 0) val.w[0] = (uint64_t)v;
    else { val.w[0] = (uint64_t)(-(int64_t)v); val = neg(val); }
  }

  void operator=(const FP192 &rhs) { val = rhs.val; }
  void operator=(uint64_t rhs) { val = U192{{rhs, 0, 0}}; }

  void setZero() { val = U192{{0, 0, 0}}; }
  void assign_no_mod(const U192 &rhs) { val = rhs; canon(val); }
  void assign_no_mod(uint64_t rhs) { val = U192{{rhs, 0, 0}}; }

  U192 value() const { return val; }

  // SPDZ layout helpers (value and MAC live in FP192x2); identity here so the
  // same templates compile, as for FP107.
  void setHigh(const FP192 &rhs) { val = rhs.val; }
  FP192 getHigh() const { return *this; }
  void setLow(const FP192 &rhs) { val = rhs.val; }
  FP192 getLow() const { return *this; }

  bool operator==(const FP192 &rhs) const { return val == rhs.val; }
  bool operator!=(const FP192 &rhs) const { return val != rhs.val; }
  bool operator==(uint64_t rhs) const { return val == U192{{rhs, 0, 0}}; }
  bool operator!=(uint64_t rhs) const { return !(*this == rhs); }

  FP192 operator+(const FP192 &rhs) const { FP192 r; r.val = add_mod(val, rhs.val); return r; }
  FP192 operator-(const FP192 &rhs) const { FP192 r; r.val = sub_mod(val, rhs.val); return r; }
  FP192 operator-() const { FP192 r; r.val = neg(val); return r; }
  FP192 operator*(const FP192 &rhs) const { FP192 r; r.val = mult_mod(val, rhs.val); return r; }

  FP192 operator+(uint64_t rhs) const { return *this + FP192(rhs); }
  FP192 operator-(uint64_t rhs) const { return *this - FP192(rhs); }
  FP192 operator*(uint64_t rhs) const { return *this * FP192(rhs); }

  // Multiply by 2^n_pos, 0 <= n_pos < 192 (COPE's bit recomposition).
  FP192 operator<<(const int n_pos) const {
    U192 t{{0, 0, 0}};
    t.w[n_pos >> 6] = 1ULL << (n_pos & 63);
    FP192 r; r.val = mult_mod(val, t); return r;
  }

  FP192 &operator+=(const FP192 &rhs) { *this = *this + rhs; return *this; }
  FP192 &operator-=(const FP192 &rhs) { *this = *this - rhs; return *this; }
  FP192 &operator*=(const FP192 &rhs) { *this = *this * rhs; return *this; }

  // No headroom above p: every add is reduced, reduce() is a no-op.
  static constexpr unsigned lazy_adds = 1;
  void add_raw(const FP192 &rhs) { val = add_mod(val, rhs.val); }
  void reduce() {}

  FP192 negate() const { return -*this; }

  // a^(p-2); returns 0 for 0.
  FP192 inv() const {
    // p - 2 = (~0, ~0 - 1, ~0 - 2)
    const U192 e{{~0ULL - 2, ~0ULL - 1, ~0ULL}};
    U192 r{{1, 0, 0}};
    for (int i = 191; i >= 0; --i) {
      r = mult_mod(r, r);
      if (e.bit(i)) r = mult_mod(r, val);
    }
    FP192 o; o.val = r; return o;
  }

  template <typename PRNG>
  void rand(PRNG &prg) {
    prg.random_data_unaligned(&val, sizeof(U192));
    canon(val);
  }

  template <typename IO> void send(IO *netio) { netio->send_data(&val, sizeof(U192)); }
  template <typename IO> void recv(IO *netio) { netio->recv_data(&val, sizeof(U192)); canon(val); }

  block hash() { return Hash::hash_for_block(&val, sizeof(U192)); }

  // The integer value of a block (< 2^128 < p): fine as a hash seed or a
  // commitment coefficient with <= 128 bits of entropy, NOT a uniform element.
  void from_block(block rhs) {
    std::memcpy(&val.w[0], &rhs, sizeof(block));
    val.w[2] = 0;
  }
  // 192 bits from two blocks (b0 whole, low half of b1): uniform up to 2^-128.
  void from_2blocks(block b0, block b1) {
    std::memcpy(&val.w[0], &b0, sizeof(block));
    std::memcpy(&val.w[2], &b1, sizeof(uint64_t));
    canon(val);
  }

  static U192 copy_compose(const U192 &basis) { return basis; }
  static uint64_t copy_compose(uint64_t basis) { return basis; }

  // ---- limb arithmetic -----------------------------------------------------

  static bool geq_p(const U192 &x) {
    // p = (~0, ~0 - 1, ~0): x >= p iff x2 = ~0 and (x1 = ~0 or (x1 = ~0-1 and x0 = ~0)).
    return x.w[2] == ~0ULL && (x.w[1] == ~0ULL || (x.w[1] == ~0ULL - 1 && x.w[0] == ~0ULL));
  }

  // The reductions below are branch-free: on random operands a carry out of
  // 2^192 happens about half the time, so branching on it mispredicts.

  // x += k * (2^64 + 1) mod 2^192; returns the carry out of 2^192.
  static uint64_t add_k(U192 &x, uint64_t k) {
    u128 t = (u128)x.w[0] + k;
    x.w[0] = (uint64_t)t;
    t = (u128)x.w[1] + k + (uint64_t)(t >> 64);
    x.w[1] = (uint64_t)t;
    t = (u128)x.w[2] + (uint64_t)(t >> 64);
    x.w[2] = (uint64_t)t;
    return (uint64_t)(t >> 64);
  }
  // x = mask ? t : x
  static void select(U192 &x, const U192 &t, uint64_t mask) {
    x.w[0] = (t.w[0] & mask) | (x.w[0] & ~mask);
    x.w[1] = (t.w[1] & mask) | (x.w[1] & ~mask);
    x.w[2] = (t.w[2] & mask) | (x.w[2] & ~mask);
  }
  // x < 2^192 -> x mod p: x >= p iff x + (2^64 + 1) carries out of 2^192.
  static void canon(U192 &x) {
    U192 t = x;
    uint64_t c = add_k(t, 1);
    select(x, t, 0 - c);
  }
  // Fold x + k*2^192 (k < 2^63) into [0, p): 2^192 = 2^64 + 1.
  static void fold(U192 &x, uint64_t k) {
    uint64_t c = add_k(x, k);
    add_k(x, c);  // x < (k + 1)*(2^64 + 1) after a carry, so no second carry
    canon(x);
  }

  static U192 add_mod(const U192 &a, const U192 &b) {
    // a + b < 2p: the result is s = a + b, or s - p = s + 2^64 + 1 (mod 2^192)
    // when a + b >= 2^192 (carry) or when s + 2^64 + 1 carries (s >= p).
    U192 s;
    u128 t = (u128)a.w[0] + b.w[0];
    s.w[0] = (uint64_t)t;
    t = (u128)a.w[1] + b.w[1] + (uint64_t)(t >> 64);
    s.w[1] = (uint64_t)t;
    t = (u128)a.w[2] + b.w[2] + (uint64_t)(t >> 64);
    s.w[2] = (uint64_t)t;
    uint64_t c = (uint64_t)(t >> 64);
    U192 r = s;
    uint64_t c2 = add_k(r, 1);
    select(s, r, 0 - (c | c2));
    return s;
  }

  static U192 sub_mod(const U192 &a, const U192 &b) {
    // d = a - b mod 2^192; on a borrow add p back, i.e. subtract 2^64 + 1.
    U192 d;
    u128 t = (u128)a.w[0] - b.w[0];
    d.w[0] = (uint64_t)t;
    t = (u128)a.w[1] - b.w[1] - (uint64_t)((t >> 64) & 1);
    d.w[1] = (uint64_t)t;
    t = (u128)a.w[2] - b.w[2] - (uint64_t)((t >> 64) & 1);
    d.w[2] = (uint64_t)t;
    uint64_t mask = 0 - (uint64_t)((t >> 64) & 1);
    uint64_t k = mask & 1;  // subtract k * (2^64 + 1)
    t = (u128)d.w[0] - k;
    d.w[0] = (uint64_t)t;
    t = (u128)d.w[1] - k - (uint64_t)((t >> 64) & 1);
    d.w[1] = (uint64_t)t;
    d.w[2] = d.w[2] - (uint64_t)((t >> 64) & 1);
    return d;
  }

  static U192 neg(const U192 &a) {
    if ((a.w[0] | a.w[1] | a.w[2]) == 0) return a;
    return sub_mod(PR, a);
  }

  static U192 mult_mod(const U192 &a, const U192 &b) {
    // 3x3 schoolbook into c[0..5].
    uint64_t c[6];
    u128 t;
    uint64_t carry;
    // row 0
    t = (u128)a.w[0] * b.w[0];                 c[0] = (uint64_t)t; carry = (uint64_t)(t >> 64);
    t = (u128)a.w[0] * b.w[1] + carry;         c[1] = (uint64_t)t; carry = (uint64_t)(t >> 64);
    t = (u128)a.w[0] * b.w[2] + carry;         c[2] = (uint64_t)t; c[3] = (uint64_t)(t >> 64);
    // row 1
    t = (u128)a.w[1] * b.w[0] + c[1];          c[1] = (uint64_t)t; carry = (uint64_t)(t >> 64);
    t = (u128)a.w[1] * b.w[1] + c[2] + carry;  c[2] = (uint64_t)t; carry = (uint64_t)(t >> 64);
    t = (u128)a.w[1] * b.w[2] + c[3] + carry;  c[3] = (uint64_t)t; c[4] = (uint64_t)(t >> 64);
    // row 2
    t = (u128)a.w[2] * b.w[0] + c[2];          c[2] = (uint64_t)t; carry = (uint64_t)(t >> 64);
    t = (u128)a.w[2] * b.w[1] + c[3] + carry;  c[3] = (uint64_t)t; carry = (uint64_t)(t >> 64);
    t = (u128)a.w[2] * b.w[2] + c[4] + carry;  c[4] = (uint64_t)t; c[5] = (uint64_t)(t >> 64);

    // Solinas: r = (c2,c1,c0) + (0,c3,c3) + (c4,c4,0) + (c5,c5,c5)  (high..low)
    U192 r;
    t = (u128)c[0] + c[3] + c[5];
    r.w[0] = (uint64_t)t;
    t = (u128)c[1] + c[3] + c[4] + c[5] + (uint64_t)(t >> 64);
    r.w[1] = (uint64_t)t;
    t = (u128)c[2] + c[4] + c[5] + (uint64_t)(t >> 64);
    r.w[2] = (uint64_t)t;
    fold(r, (uint64_t)(t >> 64));
    return r;
  }

  static std::size_t size() { return sizeof(U192); }

  // Lazy accumulator for sums of many elements (the LPN encoder, lpn.h): a
  // fourth limb absorbs the carries, one fold at the end.
  struct Acc {
    uint64_t w[4];
    explicit Acc(const FP192 &x) : w{x.val.w[0], x.val.w[1], x.val.w[2], 0} {}
    void add(const FP192 &y) {
      u128 t = (u128)w[0] + y.val.w[0];
      w[0] = (uint64_t)t;
      t = (u128)w[1] + y.val.w[1] + (uint64_t)(t >> 64);
      w[1] = (uint64_t)t;
      t = (u128)w[2] + y.val.w[2] + (uint64_t)(t >> 64);
      w[2] = (uint64_t)t;
      w[3] += (uint64_t)(t >> 64);
    }
    FP192 get() const {
      FP192 r;
      r.val = U192{{w[0], w[1], w[2]}};
      fold(r.val, w[3]);
      return r;
    }
  };
};
static_assert(sizeof(FP192) == 24, "FP192 must be 24 bytes (sent/hashed as raw bytes)");

// GGM leaves must be expanded to two blocks for this field (see spfss_*.h).
template <> struct wide_leaf<FP192> { static constexpr bool value = true; };
// Commitment-matrix coefficients come from FP::from_block on a PRP block whose
// top 32 bits also pick the row (com_matrix.h): >= 96 uniform bits per row.
template <> struct com_coef_bits<FP192> { static constexpr double value = 96.0; };

#endif  // EMP_UTILS_FIELD_FP192_H__
