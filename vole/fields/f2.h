#ifndef VOLE_FIELDS_F2_H__
#define VOLE_FIELDS_F2_H__

// Subfield VOLE with values in F_2 and keys / MACs / Delta in GF(2^128), i.e.
// Ferret's correlated OT  R = K xor b*Delta, b in {0,1}, as an instantiation
// of the generic VoleTriple<IO, FP, FPS> with
//
//   FP  = F2kKey : a GF(2^128) element whose LSB is reserved (Ferret's
//                  encoding: every sender key has LSB 0, Delta has LSB 1);
//   FPS = F2Auth : ONE block per receiver output, LSB = the value bit b and
//                  the whole block = the MAC R = K xor b*Delta.
//
// Because Delta's LSB is 1 and keys' LSBs are 0, LSB(R) = b automatically and
// XOR of blocks is the field operation on both lanes at once, so the
// accumulate / LPN / self-seeding paths are plain block XOR (16 bytes per
// output, as in emp-ot's FerretCOT). The three places where a bit value
// differs from a field value are exposed through vole_traits<F2Auth>:
//
//   unit_point_value : the MPFSS point value is the constant 1 and consumes no
//                      base pair (sender secret_sum = Delta xor sum(leaves),
//                      receiver adds nothing) -- Ferret's SPCOT;
//   mask_pairs = 128 : the consistency-check mask x* in GF(2^128) is packed
//                      from 128 base COTs, r = sum r_j X^j, M[r] = sum M_j X^j
//                      (Ferret's GaloisFieldPacking);
//   base normalisation : after COPE, sender keys get LSB 0 and receiver blocks
//                      get LSB b, which keeps R = K xor b*Delta exact because
//                      (K xor b*Delta) & ~1 = (K & ~1) xor b*(Delta & ~1).

#include "emp-tool/emp-tool.h"
#include "vole/fields/utils.h"
#include "vole/fields/fp2x128.h"

// ---- FP: GF(2^128) key with reserved LSB ----
class F2kKey : public FP2x128 {
public:
  using FP2x128::FP2x128;
  F2kKey() : FP2x128() {}
  F2kKey(const FP2x128 &v) : FP2x128(v) {}

  static const u128 lsb_mask() { return ~(u128)1; }

  // GGM leaves become keys with LSB 0 (Ferret: ggm_tree[i] & minusone).
  void from_block(block rhs) { val = from_blk(rhs) & lsb_mask(); }
  // Sender-side base keys from COPE: clear the reserved bit.
  void clear_lsb() { val &= lsb_mask(); }
};

// Delta must have LSB 1 (so that LSB(R) = b); everything else is unrestricted.
inline void restrict_delta(F2kKey &d) { d.val |= (F2kKey::u128)1; }

// ---- FPS: one block, LSB = value bit, whole block = MAC ----
class F2Auth {
public:
  using u128 = unsigned __int128;
  static constexpr uint64_t PR_bit_len       = 128;
  static constexpr uint64_t slot_stride_bits = 128;
  static constexpr uint64_t PR_num_pack      = 1;
  static constexpr uint64_t DS_byte_len      = 16;
  static constexpr unsigned lazy_adds        = 1u << 30;

  u128 val;

  F2Auth() : val(0) {}
  explicit F2Auth(u128 v) : val(v) {}

  static u128 lsb_mask() { return ~(u128)1; }

  void setZero() { val = 0; }
  // MAC = the whole block R.
  F2kKey getLow() const { F2kKey r; r.val = val; return r; }
  // value = LSB(R).
  F2kKey getHigh() const { F2kKey r; r.val = val & (u128)1; return r; }
  void setLow(const F2kKey &m) { val = (m.val & lsb_mask()) | (val & (u128)1); }
  void setHigh(const F2kKey &v) { val = (val & lsb_mask()) | (v.val & (u128)1); }
  // Base pair from COPE: (mac, value) -> (mac & ~1) | value.
  void setLowHigh(const F2kKey &m, const F2kKey &v) {
    val = (m.val & lsb_mask()) | (v.val & (u128)1);
  }
  // GGM leaf: reconstructed leaves get LSB 0 (Ferret's receiver clears too).
  void from_block(block b) { val = FP2x128::from_blk(b) & lsb_mask(); }
  void set_low_from_block(block b) { from_block(b); }

  bool operator==(const F2Auth &r) const { return val == r.val; }
  bool operator!=(const F2Auth &r) const { return val != r.val; }

  // Field operations are XOR on the whole block (both lanes at once).
  F2Auth operator+(const F2Auth &r) const { return F2Auth(val ^ r.val); }
  F2Auth operator-(const F2Auth &r) const { return F2Auth(val ^ r.val); }
  F2Auth operator+(const F2kKey &r) const { return F2Auth(val ^ r.val); }
  F2Auth operator-(const F2kKey &r) const { return F2Auth(val ^ r.val); }
  F2Auth negate() const { return *this; }
  void add_raw(const F2Auth &r) { val ^= r.val; }
  void reduce() {}
  // Broadcast multiply (only meaningful for the MAC lane; used by ComMatrixFp).
  F2Auth operator*(const F2kKey &b) const {
    return F2Auth((getLow() * b).val);
  }

  static std::size_t size() { return sizeof(u128); }
  template <typename IO> void send(IO *io) { io->send_data(&val, sizeof(u128)); }
  template <typename IO> void recv(IO *io) { io->recv_data(&val, sizeof(u128)); }
  block hash() { return Hash::hash_for_block(&val, sizeof(u128)); }
  block to_block() const { return FP2x128::to_block(val); }
};

// ---- protocol hooks for the bit instantiation ----
template <>
struct vole_traits<F2Auth> {
  static constexpr bool unit_point_value = true;
  static constexpr std::size_t mask_pairs = 128;
  using mask_t = FP2x128x2;

  // Receiver: 128 base COTs -> (r, M[r]) with r = sum r_j X^j, M[r] = sum M_j X^j.
  static mask_t pack_mask(const F2Auth *p) {
    FP2x128 r, m;
    for (std::size_t j = 0; j < 128; ++j) {
      FP2x128 xj(((unsigned __int128)1) << j, false);
      if (p[j].val & (unsigned __int128)1) r = r + xj;
      m = m + (p[j].getLow() * xj);
    }
    mask_t z; z.setLowHigh(m, r);
    return z;
  }
  // Sender: K[r] = sum K_j X^j.
  template <typename FP>
  static FP pack_mask_key(const FP *k) {
    FP2x128 s;
    for (std::size_t j = 0; j < 128; ++j) {
      FP2x128 xj(((unsigned __int128)1) << j, false);
      s = s + (k[j] * xj);
    }
    return FP(s);
  }
  // Base values are bits.
  template <typename FP, typename PRNG>
  static void rand_base_value(FP &x, PRNG &prg) {
    bool b; prg.random_bool(&b, 1);
    x.val = b ? (unsigned __int128)1 : 0;
  }
  // Sender-side base keys from COPE: reserved bit cleared.
  template <typename FP>
  static void normalize_base_keys(FP *k, std::size_t n) {
    for (std::size_t i = 0; i < n; ++i) k[i].clear_lsb();
  }
};

#endif  // VOLE_FIELDS_F2_H__
