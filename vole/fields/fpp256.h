#ifndef VOLE_FIELDS_FPP256_H__
#define VOLE_FIELDS_FPP256_H__

// Finite field F_p, p = 2^256 - 2^224 + 2^192 + 2^96 - 1 (the NIST P-256
// base field).
//
//   - Storage: four little-endian 64-bit limbs (`U256`), always reduced to
//              [0, p). No spare bit (lazy_adds = 1); the LPN encoders use Acc.
//   - Add/sub: a 257-bit sum minus p, kept iff it does not borrow; a borrow in
//              a - b adds p back. Branch-free selects.
//   - Mul:     4x4 schoolbook into eight limbs, then the FIPS 186-4 D.2.3
//              signed 32-bit word sums (T + 2S1 + 2S2 + S3 + S4 - D1 - D2 -
//              D3 - D4) carried into 64-bit limbs, one fold of the signed top
//              carry with 2^256 = 2^224 - 2^192 - 2^96 + 1, and one
//              branch-free correction by p.
//   - Inverse: Fermat, a^(p-2).
//   - Sampling: p is only ~2^-32 below 2^256, so 256 uniform bits are NOT
//              close to uniform mod p. Every uniform element is drawn from
//              384 bits and reduced (statistical distance <= p / 2^384 <
//              2^-128): rand, prg_elements, and a GGM leaf expanded to three
//              blocks (from_3blocks; leaf_blocks<FPP256> = 3).
//
// A 128-bit block maps to its integer value (< 2^128): fine as a hash seed or
// a commitment coefficient with <= 128 bits of entropy, NOT a uniform element.

#include "emp-tool/emp-tool.h"
#include "vole/fields/utils.h"
#include "vole/fields/u256.h"
#include <cstdint>
#include <cstring>

using namespace emp;

class FPP256 {
 public:
  using u128 = unsigned __int128;

  static constexpr uint64_t PR_bit_len       = 256;
  static constexpr uint64_t slot_stride_bits = 256;   // COPE decomposes Delta < 2^256
  static constexpr U256     PR               = {{~0ULL, 0x00000000FFFFFFFFULL, 0, 0xFFFFFFFF00000001ULL}};
  // Only referenced by the packed-field (PR_num_pack > 1) branches, which are
  // dead for this type; kept so those branches still compile.
  static constexpr uint64_t PR_mask     = ~0ULL;
  static constexpr uint64_t PR_num_pack = 1;
  static constexpr uint64_t DS_byte_len = 32;
  static constexpr std::size_t SAMPLE_BYTES = 48;     // bits per uniform element / 8

  U256 val;

  FPP256() : val{{0, 0, 0, 0}} {}
  FPP256(const U256 &v) : val(v) { canon(val); }   // reduces any 256-bit value
  FPP256(uint64_t v) : val{{v, 0, 0, 0}} {}
  FPP256(int v) : val{{0, 0, 0, 0}} {
    if (v >= 0) val.w[0] = (uint64_t)v;
    else { val.w[0] = (uint64_t)(-(int64_t)v); val = neg(val); }
  }

  void operator=(const FPP256 &rhs) { val = rhs.val; }
  void operator=(uint64_t rhs) { val = U256{{rhs, 0, 0, 0}}; }

  void setZero() { val = U256{{0, 0, 0, 0}}; }
  void assign_no_mod(const U256 &rhs) { val = rhs; canon(val); }
  void assign_no_mod(uint64_t rhs) { val = U256{{rhs, 0, 0, 0}}; }

  U256 value() const { return val; }

  void setHigh(const FPP256 &rhs) { val = rhs.val; }
  FPP256 getHigh() const { return *this; }
  void setLow(const FPP256 &rhs) { val = rhs.val; }
  FPP256 getLow() const { return *this; }

  bool operator==(const FPP256 &rhs) const { return val == rhs.val; }
  bool operator!=(const FPP256 &rhs) const { return val != rhs.val; }
  bool operator==(uint64_t rhs) const { return val == U256{{rhs, 0, 0, 0}}; }
  bool operator!=(uint64_t rhs) const { return !(*this == rhs); }

  FPP256 operator+(const FPP256 &rhs) const { FPP256 r; r.val = add_mod(val, rhs.val); return r; }
  FPP256 operator-(const FPP256 &rhs) const { FPP256 r; r.val = sub_mod(val, rhs.val); return r; }
  FPP256 operator-() const { FPP256 r; r.val = neg(val); return r; }
  FPP256 operator*(const FPP256 &rhs) const { FPP256 r; r.val = mult_mod(val, rhs.val); return r; }

  FPP256 operator+(uint64_t rhs) const { return *this + FPP256(rhs); }
  FPP256 operator-(uint64_t rhs) const { return *this - FPP256(rhs); }
  FPP256 operator*(uint64_t rhs) const { return *this * FPP256(rhs); }

  // Multiply by 2^n_pos, 0 <= n_pos < 256 (COPE's bit recomposition).
  FPP256 operator<<(const int n_pos) const {
    U256 t{{0, 0, 0, 0}};
    t.w[n_pos >> 6] = 1ULL << (n_pos & 63);
    FPP256 r; r.val = mult_mod(val, t); return r;
  }

  FPP256 &operator+=(const FPP256 &rhs) { *this = *this + rhs; return *this; }
  FPP256 &operator-=(const FPP256 &rhs) { *this = *this - rhs; return *this; }
  FPP256 &operator*=(const FPP256 &rhs) { *this = *this * rhs; return *this; }

  static constexpr unsigned lazy_adds = 1;
  void add_raw(const FPP256 &rhs) { val = add_mod(val, rhs.val); }
  void reduce() {}

  FPP256 negate() const { return -*this; }

  // a^(p-2); returns 0 for 0.
  FPP256 inv() const {
    const U256 e{{~0ULL - 2, 0x00000000FFFFFFFFULL, 0, 0xFFFFFFFF00000001ULL}};  // p - 2
    U256 r{{1, 0, 0, 0}};
    for (int i = 255; i >= 0; --i) {
      r = mult_mod(r, r);
      if (e.bit(i)) r = mult_mod(r, val);
    }
    FPP256 o; o.val = r; return o;
  }

  template <typename PRNG>
  void rand(PRNG &prg) {
    uint64_t b[6];
    prg.random_data_unaligned(b, sizeof(b));
    val = reduce384(b);
  }
  // n uniform elements from one PRG call (COPE's per-bit streams).
  template <typename PRNG>
  static void prg_elements(PRNG &prg, FPP256 *out, std::size_t n) {
    std::vector<uint64_t> b(6 * n);
    prg.random_data_unaligned(b.data(), b.size() * sizeof(uint64_t));
    for (std::size_t i = 0; i < n; ++i) out[i].val = reduce384(&b[6 * i]);
  }

  template <typename IO> void send(IO *netio) { netio->send_data(&val, sizeof(U256)); }
  template <typename IO> void recv(IO *netio) { netio->recv_data(&val, sizeof(U256)); canon(val); }

  block hash() { return Hash::hash_for_block(&val, sizeof(U256)); }

  void from_block(block rhs) {
    std::memcpy(&val.w[0], &rhs, sizeof(block));
    val.w[2] = val.w[3] = 0;
  }
  // 384 uniform bits from three blocks, reduced: uniform up to 2^-128.
  void from_3blocks(block b0, block b1, block b2) {
    uint64_t b[6];
    std::memcpy(&b[0], &b0, sizeof(block));
    std::memcpy(&b[2], &b1, sizeof(block));
    std::memcpy(&b[4], &b2, sizeof(block));
    val = reduce384(b);
  }

  static U256 copy_compose(const U256 &basis) { return basis; }
  static uint64_t copy_compose(uint64_t basis) { return basis; }

  // ---- limb arithmetic -----------------------------------------------------

  // x < 2^256 < 2p -> x mod p: one conditional subtract.
  static void canon(U256 &x) {
    U256 t;
    uint64_t borrow = u256::sub(t, x, PR);
    u256::select(x, t, borrow - 1);  // keep x - p iff it did not borrow
  }

  static U256 add_mod(const U256 &a, const U256 &b) {
    U256 s, t;
    uint64_t carry = u256::add(s, a, b);      // a + b < 2p < 2^257
    uint64_t borrow = u256::sub(t, s, PR);    // (carry, s) - p
    // a + b >= p iff carry or no borrow
    u256::select(s, t, 0 - (carry | (borrow ^ 1)));
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

  using i128 = __int128;
  // c = 2^256 mod p = 2^224 - 2^192 - 2^96 + 1 (< 2^224).

  // W + k * 2^256 (W < 2^256, |k| < 2^31) -> canonical residue, branch-free.
  // r = W + k*c is in (-2^255, 2^256 + 2^255), held as four limbs plus a
  // signed top carry t in {-1, 0, 1}; then one correction: + p if r < 0,
  // else - p if r >= p (one subtraction suffices: r - p < 2^256 - p < p).
  __attribute__((always_inline)) static inline U256 fold_signed(const U256 &W, int64_t k) {
    // k * c = k - k*2^96 - k*2^192 + k*2^224, as signed per-limb terms:
    // limb 0: +k, limb 1: -k*2^32, limb 3: k*2^32 - k (no multiplications).
    const i128 k32 = (i128)k << 32;
    const i128 add[4] = {(i128)k, -k32, 0, k32 - (i128)k};
    U256 r;
    i128 acc = 0;
    for (int i = 0; i < 4; ++i) {
      acc += (i128)W.w[i] + add[i];
      r.w[i] = (uint64_t)acc;
      acc >>= 64;  // arithmetic shift
    }
    const int64_t t = (int64_t)acc;
    U256 rp, rm;
    u256::add(rp, r, PR);
    const uint64_t b = u256::sub(rm, r, PR);
    const uint64_t neg = 0 - (uint64_t)(t < 0);
    const uint64_t ge = 0 - (uint64_t)((t > 0) | ((t == 0) & (b == 0)));
    u256::select(r, rm, ge & ~neg);
    u256::select(r, rp, neg);
    return r;
  }

  // Any integer < 2^512 (eight limbs) mod p: the FIPS 186-4 D.2.3 signed
  // 32-bit word sums (T + 2S1 + 2S2 + S3 + S4 - D1 - D2 - D3 - D4), carried
  // into four 64-bit limbs in 128-bit steps, leaving a signed top carry
  // k in [-5, 8] that fold_signed folds back.
  __attribute__((always_inline)) static inline U256 reduce512(const uint64_t c[8]) {
    int64_t A[16];
    for (int i = 0; i < 8; ++i) {
      A[2 * i] = (int64_t)(c[i] & 0xFFFFFFFFULL);
      A[2 * i + 1] = (int64_t)(c[i] >> 32);
    }
    int64_t s[8];
    s[0] = A[0] + A[8] + A[9] - A[11] - A[12] - A[13] - A[14];
    s[1] = A[1] + A[9] + A[10] - A[12] - A[13] - A[14] - A[15];
    s[2] = A[2] + A[10] + A[11] - A[13] - A[14] - A[15];
    s[3] = A[3] + 2 * A[11] + 2 * A[12] + A[13] - A[15] - A[8] - A[9];
    s[4] = A[4] + 2 * A[12] + 2 * A[13] + A[14] - A[9] - A[10];
    s[5] = A[5] + 2 * A[13] + 2 * A[14] + A[15] - A[10] - A[11];
    s[6] = A[6] + 3 * A[14] + 2 * A[15] + A[13] - A[8] - A[9];
    s[7] = A[7] + 3 * A[15] + A[8] - A[10] - A[11] - A[12] - A[13];
    U256 W;
    i128 acc = 0;
    for (int i = 0; i < 4; ++i) {
      acc += (i128)s[2 * i] + ((i128)s[2 * i + 1] << 32);
      W.w[i] = (uint64_t)acc;
      acc >>= 64;
    }
    return fold_signed(W, (int64_t)acc);
  }

  // 384 bits (six limbs) mod p.
  static U256 reduce384(const uint64_t b[6]) {
    uint64_t c[8] = {b[0], b[1], b[2], b[3], b[4], b[5], 0, 0};
    return reduce512(c);
  }

  static U256 mult_mod(const U256 &a, const U256 &b) {
    uint64_t c[8];
    u256::mul_wide(c, a, b);
    return reduce512(c);
  }

  static std::size_t size() { return sizeof(U256); }

  // Lazy accumulator for sums of many elements (the LPN encoders): one
  // 128-bit running sum per limb; get() forms the exact sum low + top * 2^256
  // and reduces it once (top < 2^31).
  struct Acc {
    u128 s[4];
    explicit Acc(const FPP256 &x) : s{x.val.w[0], x.val.w[1], x.val.w[2], x.val.w[3]} {}
    void add(const FPP256 &y) {
      for (int i = 0; i < 4; ++i) s[i] += y.val.w[i];
    }
    FPP256 get() const {
      U256 r;
      u128 t = s[0];
      r.w[0] = (uint64_t)t;
      for (int i = 1; i < 4; ++i) { t = s[i] + (t >> 64); r.w[i] = (uint64_t)t; }
      FPP256 o;
      o.val = fold_signed(r, (int64_t)(uint64_t)(t >> 64));
      return o;
    }
  };
};
static_assert(sizeof(FPP256) == 32, "FPP256 must be 32 bytes (sent/hashed as raw bytes)");

// GGM leaves expand to three blocks (spfss_*.h): 384 bits per element.
template <> struct wide_leaf<FPP256> { static constexpr bool value = true; };
template <> struct leaf_blocks<FPP256> { static constexpr int value = 3; };
// Commitment-matrix coefficients come from from_block on a PRP block whose top
// 32 bits also pick the row (com_matrix.h): >= 96 uniform bits per row.
template <> struct com_coef_bits<FPP256> { static constexpr double value = 96.0; };
// MPFSS check: Horner's rule (one multiplication per leaf). Unlike FP192, the
// multiplication stays out of line here, and Horner measured ~17% faster on
// the MPFSS step (same values either way).
template <> struct poly_hash_horner<FPP256> { static constexpr bool value = true; };

#endif  // VOLE_FIELDS_FPP256_H__
