#ifndef EMP_UTILS_FIELD_FP107X2_H__
#define EMP_UTILS_FIELD_FP107X2_H__

#include "emp-tool/emp-tool.h"
#include "vole/fields/utils.h"
#include "vole/fields/fp107.h"
#include <cstdint>

// FP107x2: SPDZ (value, MAC) pair over F_p, p = 2^107 - 1.
//
// Layout: two `unsigned __int128` slots — `lo` is the value share, `hi` is
// the MAC share. Total storage 32 bytes.
//
// Cannot reuse `block` (16 bytes) since 2 * 107 = 214 bits doesn't fit in a
// __m128i, so all ops are scalar (no SSE/AVX).
class FP107x2 {
public:
  using u128 = unsigned __int128;
  using i128 = __int128;

  static constexpr uint64_t PR_bit_len       = 107;
  static constexpr uint64_t slot_stride_bits = 107;             // = PR_bit_len
  static constexpr u128     PR               = (((u128)1) << 107) - 1;
  static constexpr uint64_t DS_byte_len = 32;
  static constexpr uint64_t PR_num_pack = 1;

  // Identity broadcast: large-field MAC has only one slot.
  static u128 copy_compose(u128 basis) { return basis; }
  static uint64_t copy_compose(uint64_t basis) { return basis; }

  // [0] = low slot (value), [1] = high slot (MAC).
  u128 val[2];

  FP107x2() { val[0] = 0; val[1] = 0; }

  FP107x2(const FP107x2 &inp) { val[0] = inp.val[0]; val[1] = inp.val[1]; }

  FP107x2(u128 lo, u128 hi, bool mod_it = false) {
    if(mod_it) {
      val[0] = FP107::mod(lo);
      val[1] = FP107::mod(hi);
    } else {
      val[0] = lo; val[1] = hi;
    }
  }

  FP107x2(FP107 lo, FP107 hi) {
    val[0] = lo.val; val[1] = hi.val;
  }

  FP107 operator[](bool idx) {
    FP107 r; r.val = val[idx ? 1 : 0]; return r;
  }

  FP107 getLow() const {
    FP107 r; r.val = val[0]; return r;
  }

  void setLow(FP107 lo) { val[0] = lo.val; }

  FP107 getHigh() const {
    FP107 r; r.val = val[1]; return r;
  }

  void setHigh(FP107 hi) { val[1] = hi.val; }

  void setZero() { val[0] = 0; val[1] = 0; }

  void setLowHigh(const FP107 lo, const FP107 hi) {
    val[0] = lo.val; val[1] = hi.val;
  }

  void from_block(block rhs) {
    uint64_t lo = _mm_extract_epi64(rhs, 0);
    uint64_t hi = _mm_extract_epi64(rhs, 1);
    u128 x = ((u128)hi << 64) | (u128)lo;
    val[0] = FP107::mod(x);
    val[1] = 0;
  }

  void getVal(FP107 *data) {
    data[0].val = val[0];
    data[1].val = val[1];
  }

  void operator=(const FP107x2 rhs) { val[0] = rhs.val[0]; val[1] = rhs.val[1]; }

  bool operator==(const FP107x2 rhs) const {
    return val[0] == rhs.val[0] && val[1] == rhs.val[1];
  }

  bool operator!=(const FP107x2 rhs) const {
    return !(*this == rhs);
  }

  // c = (a0 * b || a1 * b) mod PR (broadcast multiply by FP107)
  FP107x2 operator*(const FP107 b) const {
    FP107x2 res;
    res.val[0] = FP107::mult_mod(val[0], b.val);
    res.val[1] = FP107::mult_mod(val[1], b.val);
    return res;
  }

  // c = (a0 * b0 || a1 * b1) mod PR (lane-wise)
  FP107x2 operator*(const FP107x2 b) const {
    FP107x2 res;
    res.val[0] = FP107::mult_mod(val[0], b.val[0]);
    res.val[1] = FP107::mult_mod(val[1], b.val[1]);
    return res;
  }

  FP107x2 operator+(const FP107 b) const {
    FP107x2 res;
    res.val[0] = FP107::add_mod(val[0], b.val);
    res.val[1] = FP107::add_mod(val[1], b.val);
    return res;
  }

  FP107x2 operator+(const FP107x2 b) const {
    FP107x2 res;
    res.val[0] = FP107::add_mod(val[0], b.val[0]);
    res.val[1] = FP107::add_mod(val[1], b.val[1]);
    return res;
  }

  FP107x2 operator-(const FP107 b) const {
    FP107x2 res;
    res.val[0] = FP107::sub_mod(val[0], b.val);
    res.val[1] = FP107::sub_mod(val[1], b.val);
    return res;
  }

  FP107x2 operator-(const FP107x2 b) const {
    FP107x2 res;
    res.val[0] = FP107::sub_mod(val[0], b.val[0]);
    res.val[1] = FP107::sub_mod(val[1], b.val[1]);
    return res;
  }

  FP107x2 negate() const {
    FP107x2 res;
    res.val[0] = (val[0] == 0) ? 0 : (PR - val[0]);
    res.val[1] = (val[1] == 0) ? 0 : (PR - val[1]);
    return res;
  }

  static std::size_t size() { return sizeof(u128) * 2; }

  template<typename IO>
  void send(IO *netio) {
    netio->send_data(&val[0], sizeof(u128) * 2);
  }

  template<typename IO>
  void recv(IO *netio) {
    netio->recv_data(&val[0], sizeof(u128) * 2);
  }

  block hash() {
    return Hash::hash_for_block(&val[0], sizeof(u128) * 2);
  }

};

#endif // EMP_UTILS_FIELD_FP107X2_H__
