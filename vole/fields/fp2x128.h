#ifndef EMP_UTILS_FIELD_FP2X128_H__
#define EMP_UTILS_FIELD_FP2X128_H__

// Binary extension field F_{2^128} = GF(2^128), so the field-generic VOLE can
// run over "f2k" (the field the paper's evaluation uses), not just prime
// fields.
//
//   FP2x128    : a single field element (the FQ / key / MAC type).
//   FP2x128x2  : the SPDZ (value || MAC) bundle (the FPS type), two elements.
//
// Why this works through the prime-field-oriented base layer (cope.h):
//   * `val` is a `__uint128_t` so COPE's bit-decomposition of Delta
//     (delta64_to_bool: `& 1`, `>>= 1`) extracts the 128 *polynomial
//     coefficients* of Delta as true 128-bit-integer shifts.
//   * `operator<<(i)` is multiply-by-X^i in the field (not 2^i, which is 0 in
//     char 2). COPE reconstructs Delta*u = sum_i a_i * X^i, which is exactly
//     the OT-based product when Delta = sum_i Delta_i X^i.
//   * Addition is XOR; multiplication is emp's gfmul (carry-less mult + the
//     non-reflected GCM reduction, where bit i == coefficient of X^i). The
//     `__uint128_t` and the `block` fed to gfmul share byte layout on
//     little-endian, so bit i agrees across decomposition, `<<`, and gfmul.

#include "emp-tool/emp-tool.h"
#include "vole/fields/utils.h"
#include <cstdint>
#include <cstring>

using namespace emp;

class FP2x128 {
public:
  using u128 = unsigned __int128;

  // 128-bit binary field. PR_mask is all-ones (cope's `& PR_mask` is identity);
  // PR is unused for arithmetic (XOR/gfmul self-reduce) but defined for the
  // shared field interface.
  static constexpr uint64_t PR_bit_len       = 128;
  static constexpr uint64_t slot_stride_bits = 128;            // = PR_bit_len
  static constexpr u128     PR_mask          = ~(u128)0;
  static constexpr u128     PR               = ~(u128)0;
  static constexpr uint64_t PR_num_pack      = 1;
  static constexpr uint64_t DS_byte_len      = 16;

  u128 val;

  FP2x128() : val(0) {}
  FP2x128(u128 v, bool /*reduce*/ = true) : val(v) {}
  FP2x128(uint64_t v) : val((u128)v) {}
  FP2x128(int v) : val((u128)(uint64_t)v) {}

  static block to_block(u128 v) {
    block b;
    std::memcpy(&b, &v, sizeof(b));
    return b;
  }
  static u128 from_blk(block b) {
    u128 v;
    std::memcpy(&v, &b, sizeof(v));
    return v;
  }

  void operator=(const FP2x128 &rhs) { val = rhs.val; }
  void operator=(uint64_t rhs) { val = (u128)rhs; }

  void setZero() { val = 0; }
  void assign_no_mod(const u128 rhs) { val = rhs; }
  u128 value() const { return val; }

  // Single element holds either value or MAC; getHigh/getLow/setHigh/setLow are
  // identities here (the pair lives in FP2x128x2), matching FP107's surface.
  void setHigh(const u128 rhs) { val = rhs; }
  void setHigh(const FP2x128 &rhs) { val = rhs.val; }
  FP2x128 getHigh() { FP2x128 r; r.val = val; return r; }
  void setLow(const FP2x128 &rhs) { val = rhs.val; }
  FP2x128 getLow() { FP2x128 r; r.val = val; return r; }

  bool operator==(const FP2x128 &rhs) const { return val == rhs.val; }
  bool operator!=(const FP2x128 &rhs) const { return val != rhs.val; }
  bool operator==(const u128 &rhs) const { return val == rhs; }
  bool operator!=(const u128 &rhs) const { return val != rhs; }
  bool operator==(uint64_t rhs) const { return val == (u128)rhs; }
  bool operator!=(uint64_t rhs) const { return val != (u128)rhs; }

  // Addition / subtraction are XOR (characteristic 2).
  FP2x128 operator+(const FP2x128 &rhs) const { return FP2x128(val ^ rhs.val, false); }
  FP2x128 operator-(const FP2x128 &rhs) const { return FP2x128(val ^ rhs.val, false); }
  FP2x128 operator+(u128 rhs) const { return FP2x128(val ^ rhs, false); }
  FP2x128 operator-(u128 rhs) const { return FP2x128(val ^ rhs, false); }
  FP2x128 operator+(uint64_t rhs) const { return FP2x128(val ^ (u128)rhs, false); }
  FP2x128 operator-(uint64_t rhs) const { return FP2x128(val ^ (u128)rhs, false); }
  FP2x128 operator-() const { return *this; }
  FP2x128 negate() const { return *this; }
  // Addition is XOR: nothing to reduce.
  static constexpr unsigned lazy_adds = 1u << 30;
  void add_raw(const FP2x128 &rhs) { val ^= rhs.val; }
  void reduce() {}

  // Multiplication is GF(2^128) gfmul (carry-less mult + GCM reduction).
  FP2x128 operator*(const FP2x128 &rhs) const {
    block c;
    gfmul(to_block(val), to_block(rhs.val), &c);
    return FP2x128(from_blk(c), false);
  }
  FP2x128 operator*(u128 rhs) const { return *this * FP2x128(rhs, false); }
  FP2x128 operator*(uint64_t rhs) const { return *this * FP2x128((u128)rhs, false); }

  // Multiply by X^n (n in [0,127]). COPE uses this to recombine Delta*u =
  // sum_n a_n * X^n. X^n for n<128 is the element with bit n set.
  FP2x128 operator<<(const int n_pos) const {
    return *this * FP2x128(((u128)1) << n_pos, false);
  }

  // Multiplicative inverse a^(2^128 - 2) = prod_{i=1}^{127} a^(2^i) (0 -> 0).
  FP2x128 inv() const {
    FP2x128 t = *this, r(1);
    for (int i = 1; i < 128; ++i) { t = t * t; r = r * t; }
    return r;
  }

  FP2x128 &operator+=(const FP2x128 &rhs) { val ^= rhs.val; return *this; }
  FP2x128 &operator-=(const FP2x128 &rhs) { val ^= rhs.val; return *this; }
  FP2x128 &operator*=(const FP2x128 &rhs) { *this = *this * rhs; return *this; }

  template <typename PRNG>
  void rand(PRNG &prg) {
    prg.random_data(&val, sizeof(u128));  // every 128-bit string is a valid elt
  }

  template <typename IO>
  void send(IO *netio) { netio->send_data(&val, sizeof(u128)); }
  template <typename IO>
  void recv(IO *netio) { netio->recv_data(&val, sizeof(u128)); }

  block hash() { return Hash::hash_for_block(&val, sizeof(u128)); }

  void from_block(block rhs) { val = from_blk(rhs); }

  static u128 copy_compose(u128 basis) { return basis; }
  static uint64_t copy_compose(uint64_t basis) { return basis; }


  static std::size_t size() { return sizeof(u128); }
};

// SPDZ (value || MAC) pair over GF(2^128). Two 128-bit slots (32 bytes); a
// single block (16 bytes) cannot hold both, so ops are scalar.
//   val[0] = low slot (MAC), val[1] = high slot (value).
// (setLowHigh(lo, hi) is called by the base sVOLE as setLowHigh(mac, value),
//  so getLow()=MAC, getHigh()=value, matching the rest of the VOLE code.)
class FP2x128x2 {
public:
  using u128 = unsigned __int128;

  static constexpr uint64_t PR_bit_len       = 128;
  static constexpr uint64_t slot_stride_bits = 128;
  static constexpr u128     PR               = ~(u128)0;
  static constexpr uint64_t PR_num_pack      = 1;
  static constexpr uint64_t DS_byte_len      = 32;

  static u128 copy_compose(u128 basis) { return basis; }
  static uint64_t copy_compose(uint64_t basis) { return basis; }

  u128 val[2];

  FP2x128x2() { val[0] = 0; val[1] = 0; }
  FP2x128x2(const FP2x128x2 &inp) { val[0] = inp.val[0]; val[1] = inp.val[1]; }
  FP2x128x2(u128 lo, u128 hi, bool /*mod_it*/ = false) { val[0] = lo; val[1] = hi; }
  FP2x128x2(FP2x128 lo, FP2x128 hi) { val[0] = lo.val; val[1] = hi.val; }

  FP2x128 operator[](bool idx) { FP2x128 r; r.val = val[idx ? 1 : 0]; return r; }

  FP2x128 getLow() const { FP2x128 r; r.val = val[0]; return r; }
  void setLow(FP2x128 lo) { val[0] = lo.val; }
  FP2x128 getHigh() const { FP2x128 r; r.val = val[1]; return r; }
  void setHigh(FP2x128 hi) { val[1] = hi.val; }
  void setZero() { val[0] = 0; val[1] = 0; }
  void setLowHigh(const FP2x128 lo, const FP2x128 hi) { val[0] = lo.val; val[1] = hi.val; }

  void from_block(block rhs) {
    val[0] = FP2x128::from_blk(rhs);  // low (MAC) slot
    val[1] = 0;
  }

  void getVal(FP2x128 *data) { data[0].val = val[0]; data[1].val = val[1]; }

  void operator=(const FP2x128x2 rhs) { val[0] = rhs.val[0]; val[1] = rhs.val[1]; }

  bool operator==(const FP2x128x2 rhs) const {
    return val[0] == rhs.val[0] && val[1] == rhs.val[1];
  }
  bool operator!=(const FP2x128x2 rhs) const { return !(*this == rhs); }

  // Broadcast multiply by a single element (gfmul each slot).
  FP2x128x2 operator*(const FP2x128 b) const {
    FP2x128x2 res;
    res.val[0] = (getLow() * b).val;
    res.val[1] = (getHigh() * b).val;
    return res;
  }
  // Lane-wise multiply.
  FP2x128x2 operator*(const FP2x128x2 b) const {
    FP2x128x2 res;
    res.val[0] = (getLow() * b.getLow()).val;
    res.val[1] = (getHigh() * b.getHigh()).val;
    return res;
  }

  // Addition is XOR. operator+(FP) XORs into BOTH slots (the value slot is
  // overwritten by set_vec_x afterwards, matching FP61x2 / FP107x2).
  FP2x128x2 operator+(const FP2x128 b) const {
    return FP2x128x2(val[0] ^ b.val, val[1] ^ b.val);
  }
  FP2x128x2 operator+(const FP2x128x2 b) const {
    return FP2x128x2(val[0] ^ b.val[0], val[1] ^ b.val[1]);
  }
  FP2x128x2 operator-(const FP2x128 b) const { return *this + b; }
  FP2x128x2 operator-(const FP2x128x2 b) const { return *this + b; }
  FP2x128x2 negate() const { return *this; }
  static constexpr unsigned lazy_adds = 1u << 30;
  void add_raw(const FP2x128x2 &rhs) { val[0] ^= rhs.val[0]; val[1] ^= rhs.val[1]; }
  void reduce() {}
  void set_low_from_block(block b) { val[0] = FP2x128::from_blk(b); val[1] = 0; }

  static std::size_t size() { return sizeof(u128) * 2; }

  template <typename IO>
  void send(IO *netio) { netio->send_data(&val[0], sizeof(u128) * 2); }
  template <typename IO>
  void recv(IO *netio) { netio->recv_data(&val[0], sizeof(u128) * 2); }

  block hash() { return Hash::hash_for_block(&val[0], sizeof(u128) * 2); }

};

#endif  // EMP_UTILS_FIELD_FP2X128_H__
