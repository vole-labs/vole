#ifndef VOLE_FIELDS_FIELD_PAIR_H__
#define VOLE_FIELDS_FIELD_PAIR_H__

// FieldPair<F>: SPDZ (value, MAC) pair of two F slots for the limb-storage
// fields (FP128, FP25519, FPP256), the same layout and operations as FP192x2:
// val[0] is the low slot (the MAC in a VOLE bundle), val[1] the high slot
// (the value); all operations are lane-wise.

#include "emp-tool/emp-tool.h"
#include "vole/fields/utils.h"
#include <cstdint>

using namespace emp;

template <typename F>
class FieldPair {
public:
  using limb_t = decltype(F().val);
  static constexpr uint64_t PR_bit_len       = F::PR_bit_len;
  static constexpr uint64_t slot_stride_bits = F::slot_stride_bits;
  static constexpr uint64_t DS_byte_len = 2 * sizeof(F);
  static constexpr uint64_t PR_num_pack = 1;

  static limb_t copy_compose(const limb_t &basis) { return basis; }
  static uint64_t copy_compose(uint64_t basis) { return basis; }

  F val[2];

  FieldPair() {}
  FieldPair(const FieldPair &inp) { val[0] = inp.val[0]; val[1] = inp.val[1]; }
  FieldPair(F lo, F hi) { val[0] = lo; val[1] = hi; }

  F operator[](bool idx) { return val[idx ? 1 : 0]; }

  F getLow() const { return val[0]; }
  void setLow(F lo) { val[0] = lo; }
  F getHigh() const { return val[1]; }
  void setHigh(F hi) { val[1] = hi; }
  void setZero() { val[0].setZero(); val[1].setZero(); }
  void setLowHigh(const F lo, const F hi) { val[0] = lo; val[1] = hi; }

  void from_block(block rhs) { val[0].from_block(rhs); val[1].setZero(); }
  void set_low_from_block(block b) { val[0].from_block(b); val[1].setZero(); }
  void set_low_from_2blocks(block b0, block b1)
    requires requires(F f) { f.from_2blocks(b0, b1); } {
    val[0].from_2blocks(b0, b1); val[1].setZero();
  }
  void set_low_from_3blocks(block b0, block b1, block b2)
    requires requires(F f) { f.from_3blocks(b0, b1, b2); } {
    val[0].from_3blocks(b0, b1, b2); val[1].setZero();
  }

  void getVal(F *data) { data[0] = val[0]; data[1] = val[1]; }

  void operator=(const FieldPair rhs) { val[0] = rhs.val[0]; val[1] = rhs.val[1]; }

  bool operator==(const FieldPair rhs) const { return val[0] == rhs.val[0] && val[1] == rhs.val[1]; }
  bool operator!=(const FieldPair rhs) const { return !(*this == rhs); }

  FieldPair operator*(const F b) const { return FieldPair(val[0] * b, val[1] * b); }
  FieldPair operator*(const FieldPair b) const { return FieldPair(val[0] * b.val[0], val[1] * b.val[1]); }
  FieldPair operator+(const F b) const { return FieldPair(val[0] + b, val[1] + b); }
  FieldPair operator+(const FieldPair b) const { return FieldPair(val[0] + b.val[0], val[1] + b.val[1]); }
  FieldPair operator-(const F b) const { return FieldPair(val[0] - b, val[1] - b); }
  FieldPair operator-(const FieldPair b) const { return FieldPair(val[0] - b.val[0], val[1] - b.val[1]); }

  static constexpr unsigned lazy_adds = 1;
  void add_raw(const FieldPair &rhs) { val[0] = val[0] + rhs.val[0]; val[1] = val[1] + rhs.val[1]; }
  void reduce() {}

  FieldPair negate() const { return FieldPair(val[0].negate(), val[1].negate()); }

  static std::size_t size() { return 2 * sizeof(F); }

  template <typename IO> void send(IO *netio) { netio->send_data(&val[0], 2 * sizeof(F)); }
  template <typename IO> void recv(IO *netio) {
    netio->recv_data(&val[0], 2 * sizeof(F));
    F::canon(val[0].val);
    F::canon(val[1].val);
  }

  block hash() { return Hash::hash_for_block(&val[0], 2 * sizeof(F)); }

  // Lane-wise lazy accumulator (see F::Acc).
  struct Acc {
    typename F::Acc a0, a1;
    explicit Acc(const FieldPair &x) : a0(x.val[0]), a1(x.val[1]) {}
    void add(const FieldPair &y) { a0.add(y.val[0]); a1.add(y.val[1]); }
    FieldPair get() const { return FieldPair(a0.get(), a1.get()); }
  };
};

#endif  // VOLE_FIELDS_FIELD_PAIR_H__
