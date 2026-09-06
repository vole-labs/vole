#ifndef EMP_UTILS_FIELD_Z2K_H__
#define EMP_UTILS_FIELD_Z2K_H__

// VOLE over the ring Z_{2^k} (SPDZ2k-style), so the field-generic committed
// VOLE can run over a ring, not only over prime / GF(2^128) fields.
//
//   Z2k64    : a single ring element (the key / MAC type, the FP role).
//   Z2k64x2  : the SPDZ (value || MAC) bundle (the FPS role), two elements.
//
// Parameters (SPDZ2k-style): the authenticated value lives in
// Z_{2^k} with k = 64, the MAC is computed in Z_{2^{k+s}} with k+s = 128, and
// Delta lives in Z_{2^s} with s = 64. Concretely:
//   * `val` is an unsigned __int128, so ALL arithmetic (+, -, *) is the native
//     two's-complement operation, i.e. reduction mod 2^128 is the hardware
//     overflow -- there is no explicit modular reduction (unlike the Mersenne
//     fields) and no carry-less multiply (unlike GF(2^128)).
//   * `operator<<(i)` is multiply-by-2^i (native shift). COPE reconstructs
//     Delta*u = sum_i a_i * 2^i, exactly the OT-based product when
//     Delta = sum_i Delta_i 2^i. This is the Gilboa / bit-decomposition product
//     that already underlies cope.h.
//   * `slot_stride_bits = 64`, so COPE decomposes only the low 64 bits of
//     Delta -- Delta MUST be confined to Z_{2^64}. A freshly sampled Delta is
//     clamped by restrict_delta() below (the s=64 MAC-soundness margin is the
//     high 64 bits of the 128-bit MAC ring).
//
// Why this is correct over a ring: the silent VOLE (GGM/SPFSS expand, COPE base
// sVOLE, accumulate + sparse LPN) is a linear map plus a single scalar
// multiply-by-Delta, all of which are defined over any commutative ring. No step
// needs a field inverse. The malicious consistency checks (base sVOLE check,
// MPFSS batch check) are universal-hash checks whose soundness over Z_{2^k}
// rests on the extra s bits (the usual SPDZ2k heuristic).

#include "emp-tool/emp-tool.h"
#include "vole/fields/utils.h"
#include <cstdint>
#include <cstring>

using namespace emp;

class Z2k64 {
public:
  using u128 = unsigned __int128;

  // k + s = 128 ring; value in Z_{2^64}, Delta in Z_{2^64}, MAC in Z_{2^128}.
  static constexpr uint64_t PR_bit_len       = 128;
  // COPE decomposes this many low bits of Delta (the s = 64 Delta sub-ring).
  static constexpr uint64_t slot_stride_bits = 64;
  // `& PR_mask` must be identity inside COPE (no field reduction on the ring).
  static constexpr u128     PR_mask          = ~(u128)0;
  static constexpr u128     PR               = ~(u128)0;
  static constexpr uint64_t PR_num_pack      = 1;
  static constexpr uint64_t DS_byte_len      = 16;
  // Admissible Delta space: low slot_stride_bits bits.
  static constexpr u128     DELTA_MASK       = ((u128)1 << slot_stride_bits) - 1;

  u128 val;

  Z2k64() : val(0) {}
  Z2k64(u128 v, bool /*reduce*/ = true) : val(v) {}
  Z2k64(uint64_t v) : val((u128)v) {}
  Z2k64(int v) : val((u128)(uint64_t)v) {}

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

  void operator=(const Z2k64 &rhs) { val = rhs.val; }
  void operator=(uint64_t rhs) { val = (u128)rhs; }

  void setZero() { val = 0; }
  void assign_no_mod(const u128 rhs) { val = rhs; }
  u128 value() const { return val; }

  // A single element holds either a value or a MAC; getHigh/getLow/setHigh/
  // setLow are identities here (the pair lives in Z2k64x2), matching the FP role
  // of FP2x128 / FP107.
  void setHigh(const u128 rhs) { val = rhs; }
  void setHigh(const Z2k64 &rhs) { val = rhs.val; }
  Z2k64 getHigh() { Z2k64 r; r.val = val; return r; }
  void setLow(const Z2k64 &rhs) { val = rhs.val; }
  Z2k64 getLow() { Z2k64 r; r.val = val; return r; }

  bool operator==(const Z2k64 &rhs) const { return val == rhs.val; }
  bool operator!=(const Z2k64 &rhs) const { return val != rhs.val; }
  bool operator==(const u128 &rhs) const { return val == rhs; }
  bool operator!=(const u128 &rhs) const { return val != rhs; }
  bool operator==(uint64_t rhs) const { return val == (u128)rhs; }
  bool operator!=(uint64_t rhs) const { return val != (u128)rhs; }

  // Addition / subtraction are native ring ops (mod 2^128 via overflow).
  Z2k64 operator+(const Z2k64 &rhs) const { return Z2k64(val + rhs.val, false); }
  Z2k64 operator-(const Z2k64 &rhs) const { return Z2k64(val - rhs.val, false); }
  Z2k64 operator+(u128 rhs) const { return Z2k64(val + rhs, false); }
  Z2k64 operator-(u128 rhs) const { return Z2k64(val - rhs, false); }
  Z2k64 operator+(uint64_t rhs) const { return Z2k64(val + (u128)rhs, false); }
  Z2k64 operator-(uint64_t rhs) const { return Z2k64(val - (u128)rhs, false); }
  Z2k64 operator-() const { return Z2k64((u128)0 - val, false); }
  Z2k64 negate() const { return Z2k64((u128)0 - val, false); }

  // Multiplication is native ring multiply (low 128 bits of the product).
  Z2k64 operator*(const Z2k64 &rhs) const { return Z2k64(val * rhs.val, false); }
  Z2k64 operator*(u128 rhs) const { return Z2k64(val * rhs, false); }
  Z2k64 operator*(uint64_t rhs) const { return Z2k64(val * (u128)rhs, false); }

  // Multiply by 2^n (n in [0,127]). COPE uses this to recombine
  // Delta*u = sum_n a_n * 2^n.
  Z2k64 operator<<(const int n_pos) const { return Z2k64(val << n_pos, false); }

  Z2k64 &operator+=(const Z2k64 &rhs) { val += rhs.val; return *this; }
  Z2k64 &operator-=(const Z2k64 &rhs) { val -= rhs.val; return *this; }
  Z2k64 &operator*=(const Z2k64 &rhs) { val *= rhs.val; return *this; }

  template <typename PRNG>
  void rand(PRNG &prg) {
    prg.random_data(&val, sizeof(u128));  // value side is full-width random
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

// SPDZ (value || MAC) pair over the ring. Two 128-bit slots (32 bytes).
//   val[0] = low slot (MAC), val[1] = high slot (value).
// (setLowHigh(lo, hi) is called by the base sVOLE as setLowHigh(mac, value),
//  so getLow()=MAC, getHigh()=value, matching the rest of the VOLE code.)
class Z2k64x2 {
public:
  using u128 = unsigned __int128;

  static constexpr uint64_t PR_bit_len       = 128;
  static constexpr uint64_t slot_stride_bits = 64;
  static constexpr u128     PR               = ~(u128)0;
  static constexpr uint64_t PR_num_pack      = 1;
  static constexpr uint64_t DS_byte_len      = 32;

  static u128 copy_compose(u128 basis) { return basis; }
  static uint64_t copy_compose(uint64_t basis) { return basis; }

  u128 val[2];

  Z2k64x2() { val[0] = 0; val[1] = 0; }
  Z2k64x2(const Z2k64x2 &inp) { val[0] = inp.val[0]; val[1] = inp.val[1]; }
  Z2k64x2(u128 lo, u128 hi, bool /*mod_it*/ = false) { val[0] = lo; val[1] = hi; }
  Z2k64x2(Z2k64 lo, Z2k64 hi) { val[0] = lo.val; val[1] = hi.val; }

  Z2k64 operator[](bool idx) { Z2k64 r; r.val = val[idx ? 1 : 0]; return r; }

  Z2k64 getLow() const { Z2k64 r; r.val = val[0]; return r; }
  void setLow(Z2k64 lo) { val[0] = lo.val; }
  Z2k64 getHigh() const { Z2k64 r; r.val = val[1]; return r; }
  void setHigh(Z2k64 hi) { val[1] = hi.val; }
  void setZero() { val[0] = 0; val[1] = 0; }
  void setLowHigh(const Z2k64 lo, const Z2k64 hi) { val[0] = lo.val; val[1] = hi.val; }

  void from_block(block rhs) {
    val[0] = Z2k64::from_blk(rhs);  // low (MAC) slot
    val[1] = 0;
  }

  void getVal(Z2k64 *data) { data[0].val = val[0]; data[1].val = val[1]; }

  void operator=(const Z2k64x2 rhs) { val[0] = rhs.val[0]; val[1] = rhs.val[1]; }

  bool operator==(const Z2k64x2 rhs) const {
    return val[0] == rhs.val[0] && val[1] == rhs.val[1];
  }
  bool operator!=(const Z2k64x2 rhs) const { return !(*this == rhs); }

  // Broadcast multiply by a single element (each slot independently).
  Z2k64x2 operator*(const Z2k64 b) const {
    return Z2k64x2(val[0] * b.val, val[1] * b.val);
  }
  // Lane-wise multiply.
  Z2k64x2 operator*(const Z2k64x2 b) const {
    return Z2k64x2(val[0] * b.val[0], val[1] * b.val[1]);
  }

  // operator+(FP) adds into BOTH slots (the value slot is overwritten by
  // set_vec_x afterwards, matching FP61x2 / FP2x128x2).
  Z2k64x2 operator+(const Z2k64 b) const {
    return Z2k64x2(val[0] + b.val, val[1] + b.val);
  }
  Z2k64x2 operator+(const Z2k64x2 b) const {
    return Z2k64x2(val[0] + b.val[0], val[1] + b.val[1]);
  }
  Z2k64x2 operator-(const Z2k64 b) const {
    return Z2k64x2(val[0] - b.val, val[1] - b.val);
  }
  Z2k64x2 operator-(const Z2k64x2 b) const {
    return Z2k64x2(val[0] - b.val[0], val[1] - b.val[1]);
  }
  Z2k64x2 negate() const { return Z2k64x2((u128)0 - val[0], (u128)0 - val[1]); }

  static std::size_t size() { return sizeof(u128) * 2; }

  template <typename IO>
  void send(IO *netio) { netio->send_data(&val[0], sizeof(u128) * 2); }
  template <typename IO>
  void recv(IO *netio) { netio->recv_data(&val[0], sizeof(u128) * 2); }

  block hash() { return Hash::hash_for_block(&val[0], sizeof(u128) * 2); }

};

// Clamp a freshly sampled Delta to the authentication sub-ring Z_{2^s}. COPE
// decomposes only the low slot_stride_bits bits of Delta, so Delta must not have
// higher bits set (otherwise the verifier's checks use a Delta the MAC was never
// formed with). No-op for the field types (see the generic in utils.h).
inline void restrict_delta(Z2k64 &d) { d.val &= Z2k64::DELTA_MASK; }

// Select the ring (MozZ2karella) malicious checks for the Z_{2^k} key type: a
// GGM-tree tag check over GF(2^128) plus a binary-coefficient linear-combination
// check, instead of the field random-coefficient check. See vole/fields/utils.h.
template<> struct field_traits<Z2k64> { static constexpr bool is_ring = true; };

#endif  // EMP_UTILS_FIELD_Z2K_H__
