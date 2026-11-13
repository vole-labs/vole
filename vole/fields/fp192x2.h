#ifndef EMP_UTILS_FIELD_FP192X2_H__
#define EMP_UTILS_FIELD_FP192X2_H__

#include "emp-tool/emp-tool.h"
#include "vole/fields/utils.h"
#include "vole/fields/fp192.h"
#include <cstdint>

// FP192x2: SPDZ (value, MAC) pair over F_p, p = 2^192 - 2^64 - 1.
//
// Layout: two FP192 slots -- val[0] is the value share, val[1] the MAC share.
// Total storage 48 bytes; all ops are scalar limb arithmetic.
class FP192x2 {
public:
  static constexpr uint64_t PR_bit_len       = 192;
  static constexpr uint64_t slot_stride_bits = 192;
  static constexpr uint64_t DS_byte_len = 48;
  static constexpr uint64_t PR_num_pack = 1;

  static U192 copy_compose(const U192 &basis) { return basis; }
  static uint64_t copy_compose(uint64_t basis) { return basis; }

  // [0] = low slot (value), [1] = high slot (MAC).
  FP192 val[2];

  FP192x2() {}
  FP192x2(const FP192x2 &inp) { val[0] = inp.val[0]; val[1] = inp.val[1]; }
  FP192x2(FP192 lo, FP192 hi) { val[0] = lo; val[1] = hi; }

  FP192 operator[](bool idx) { return val[idx ? 1 : 0]; }

  FP192 getLow() const { return val[0]; }
  void setLow(FP192 lo) { val[0] = lo; }
  FP192 getHigh() const { return val[1]; }
  void setHigh(FP192 hi) { val[1] = hi; }
  void setZero() { val[0].setZero(); val[1].setZero(); }
  void setLowHigh(const FP192 lo, const FP192 hi) { val[0] = lo; val[1] = hi; }

  void from_block(block rhs) { val[0].from_block(rhs); val[1].setZero(); }
  void set_low_from_block(block b) { val[0].from_block(b); val[1].setZero(); }
  void set_low_from_2blocks(block b0, block b1) { val[0].from_2blocks(b0, b1); val[1].setZero(); }

  void getVal(FP192 *data) { data[0] = val[0]; data[1] = val[1]; }

  void operator=(const FP192x2 rhs) { val[0] = rhs.val[0]; val[1] = rhs.val[1]; }

  bool operator==(const FP192x2 rhs) const { return val[0] == rhs.val[0] && val[1] == rhs.val[1]; }
  bool operator!=(const FP192x2 rhs) const { return !(*this == rhs); }

  // broadcast multiply by a scalar
  FP192x2 operator*(const FP192 b) const { return FP192x2(val[0] * b, val[1] * b); }
  // lane-wise
  FP192x2 operator*(const FP192x2 b) const { return FP192x2(val[0] * b.val[0], val[1] * b.val[1]); }
  FP192x2 operator+(const FP192 b) const { return FP192x2(val[0] + b, val[1] + b); }
  FP192x2 operator+(const FP192x2 b) const { return FP192x2(val[0] + b.val[0], val[1] + b.val[1]); }
  FP192x2 operator-(const FP192 b) const { return FP192x2(val[0] - b, val[1] - b); }
  FP192x2 operator-(const FP192x2 b) const { return FP192x2(val[0] - b.val[0], val[1] - b.val[1]); }

  static constexpr unsigned lazy_adds = 1;
  void add_raw(const FP192x2 &rhs) { val[0] = val[0] + rhs.val[0]; val[1] = val[1] + rhs.val[1]; }
  void reduce() {}

  FP192x2 negate() const { return FP192x2(val[0].negate(), val[1].negate()); }

  static std::size_t size() { return 2 * sizeof(FP192); }

  template <typename IO> void send(IO *netio) { netio->send_data(&val[0], 2 * sizeof(FP192)); }
  template <typename IO> void recv(IO *netio) {
    netio->recv_data(&val[0], 2 * sizeof(FP192));
    FP192::canon(val[0].val);
    FP192::canon(val[1].val);
  }

  block hash() { return Hash::hash_for_block(&val[0], 2 * sizeof(FP192)); }

  // Lane-wise lazy accumulator (see FP192::Acc).
  struct Acc {
    FP192::Acc a0, a1;
    explicit Acc(const FP192x2 &x) : a0(x.val[0]), a1(x.val[1]) {}
    void add(const FP192x2 &y) { a0.add(y.val[0]); a1.add(y.val[1]); }
    FP192x2 get() const { return FP192x2(a0.get(), a1.get()); }
  };
};
static_assert(sizeof(FP192x2) == 48, "FP192x2 must be 48 bytes");

#endif // EMP_UTILS_FIELD_FP192X2_H__
