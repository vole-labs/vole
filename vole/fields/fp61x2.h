#ifndef EMP_UTILS_FIELD_FP61_X2_H__
#define EMP_UTILS_FIELD_FP61_X2_H__

#include "emp-tool/emp-tool.h"
#include "vole/fields/utils.h"
#include "vole/fields/fp61.h"
#include <cstdint>

// define finite field Fp, p = 2^61 - 1
// define FP61x2 as Fp||Fp
class FP61x2 {
public:
  static constexpr uint64_t PR_bit_len       = 61;
  static constexpr uint64_t slot_stride_bits = 61;   // = PR_bit_len
  static constexpr uint64_t PR               = 2305843009213693951;  // 2^61 - 1
  static constexpr uint64_t DS_byte_len      = 16;
  // Single-slot pack: getLow / getHigh are the two roles, but the SIMD-MAC
  // amplification factor is 1 (each MAC slot is the full 61-bit value, no
  // packing). Exposed so field-generic templates can iterate `< PR_num_pack`.
  static constexpr uint64_t PR_num_pack      = 1ULL;
  // block constant via accessor (see FP61::DoublePR) to keep the header
  // single-definition across translation units (no file-scope static -> no ODR
  // clash when multiple .cpp files include this header).
  static block DoublePR() {
    return makeBlock(2305843009213693951, 2305843009213693951);
  }

  // Identity broadcast: large-field MAC has only one slot.
  static uint64_t copy_compose(uint64_t basis) { return basis; }

  block val;

  FP61x2() { val = zero_block; }

  FP61x2(const FP61x2 &inp) : val(inp.val) {}

  FP61x2(block val_, bool mod_it = false) : val(val_) {
    if(mod_it == true) {
      uint64_t low = _mm_extract_epi64(val, 0);
      uint64_t high = _mm_extract_epi64(val, 1);
      low = mod(low);
      high = mod(high);
      val = makeBlock(high, low);
    }
  }

  FP61x2(uint64_t low, uint64_t high, bool mod_it = false) {
    if(mod_it == true) {
      low = mod(low);
      high = mod(high);
    }
    val = makeBlock(high, low);
  }

  FP61x2(FP61 low, FP61 high) {
    val = makeBlock(high.val, low.val);
  }

  FP61 operator[](bool idx) {
    return FP61((idx == true)?_mm_extract_epi64(this->val, 1):
          _mm_extract_epi64(this->val, 0));
  }

  FP61 getLow() const {
    return FP61(_mm_extract_epi64(this->val, 0), false);
  }

  void setLow(FP61 low) {
    uint64_t high = _mm_extract_epi64(this->val, 1);
    this->val = makeBlock(high, low.val);
  }

  FP61 getHigh() const {
    return FP61(_mm_extract_epi64(this->val, 1), false);
  }

  void setHigh(FP61 high) {
    uint64_t low = _mm_extract_epi64(this->val, 0);
    val = makeBlock(high.val, low);
  }

  void setZero() {
    val = zero_block;
  }

  void setLowHigh(const FP61 low, const FP61 high) {
    val = makeBlock(high.val, low.val);
  }

  void from_block(block rhs) {
    val = makeBlock(mod(_mm_extract_epi64(rhs, 1)), mod(_mm_extract_epi64(rhs, 0)));
  }

  void getVal(FP61 *data) {
    data[0] = FP61(_mm_extract_epi64(val, 0));
    data[1] = FP61(_mm_extract_epi64(val, 1));
  }

  void operator=(const block rhs) { this->val = rhs; }

  void operator=(const FP61x2 rhs) { this->val = rhs.val; }

  bool operator==(const block rhs) const {
    __m128i vcmp = _mm_xor_si128(this->val, rhs);
    if(_mm_testz_si128(vcmp, vcmp))
      return true;
    return false;
  }

  bool operator==(const FP61x2 rhs) const {
    __m128i vcmp = _mm_xor_si128(this->val, rhs.val);
    if(_mm_testz_si128(vcmp, vcmp))
      return true;
    return false;
  }

  bool operator!=(const block rhs) const {
    __m128i vcmp = _mm_xor_si128(this->val, rhs);
    if(_mm_testz_si128(vcmp, vcmp))
      return false;
    return true;
  }

  bool operator!=(const FP61x2 rhs) const {
    __m128i vcmp = _mm_xor_si128(this->val, rhs.val);
    if(_mm_testz_si128(vcmp, vcmp))
      return false;
    return true;
  }

  uint64_t mod(uint64_t x) {
    uint64_t i = (x & PR) + (x >> PR_bit_len);
    return (i >= PR) ? i - PR : i;
  }

  // c = (a1 * b || a2 * b) mod PR
  FP61x2 operator*(const FP61 b) const {
    FP61x2 res(*this);
    uint64_t H = _mm_extract_epi64(this->val, 1);
    uint64_t L = _mm_extract_epi64(this->val, 0);
    block bs[2];
    uint64_t *is = (uint64_t *)(bs);
    is[1] = mul64(H, b.val, (uint64_t *)(is + 3));
    is[0] = mul64(L, b.val, (uint64_t *)(is + 2));
    block t1 = bs[0] & DoublePR();
    block t2 = _mm_srli_epi64(bs[0], PR_bit_len) ^
               _mm_slli_epi64(bs[1], 64 - PR_bit_len);
    res.val = _mm_add_epi64(t1, t2);
    res.val = vec_partial_mod(res.val);
    return res;
  }

  // c = (a1 * b || a2 * b) mod PR
  FP61x2 operator*(const FP61x2 b) const {
    FP61x2 res(*this);
    uint64_t H = _mm_extract_epi64(this->val, 1);
    uint64_t L = _mm_extract_epi64(this->val, 0);
    block bs[2];
    uint64_t *is = (uint64_t *)(bs);
    is[1] = mul64(H, b.getHigh().val, (uint64_t *)(is + 3));
    is[0] = mul64(L, b.getLow().val, (uint64_t *)(is + 2));
    block t1 = bs[0] & DoublePR();
    block t2 = _mm_srli_epi64(bs[0], PR_bit_len) ^
               _mm_slli_epi64(bs[1], 64 - PR_bit_len);
    res.val = _mm_add_epi64(t1, t2);
    res.val = vec_partial_mod(res.val);
    return res;
  }

  // c = (a1 + b || a2 + b) mod PR
  FP61x2 operator+(const FP61 b) const {
    FP61x2 res(*this);
    res.val = _mm_add_epi64(
        res.val, _mm_set_epi64((__m64)(b.val), (__m64)(b.val)));
    res.val = vec_partial_mod(res.val);
    return res;
  }

  // c = (a1 + b1 || a2 + b2) mod PR
  FP61x2 operator+(const FP61x2 b) const {
    FP61x2 res(*this);
    res.val = _mm_add_epi64(res.val, b.val);
    res.val = vec_partial_mod(res.val);
    return res;
  }

  // c = (a1 - b || a2 - b) mod PR
  FP61x2 operator-(const FP61 b) const {
    FP61x2 res(*this);
    __m64 b_neg = (__m64)(PR - b.val);
    res.val = _mm_add_epi64(
        res.val, _mm_set_epi64(b_neg, b_neg));
    res.val = vec_partial_mod(res.val);
    return res;
  }

  // c = (a1 - b1 || a2 - b2) mod PR
  FP61x2 operator-(const FP61x2 b) const {
    FP61x2 res(*this);
    block c = _mm_sub_epi64(DoublePR(), b.val);
    res.val = _mm_add_epi64(res.val, c);
    res.val = vec_partial_mod(res.val);
    return res;
  }

  FP61x2 negate() const {
    FP61x2 res(*this);
    uint64_t low = _mm_extract_epi64(res.val, 0);
    uint64_t high = _mm_extract_epi64(res.val, 1);
    if(low != 0) low = PR - low;
    if(high != 0) high = PR - high;
    res.val = makeBlock(high, low);
    return res;
  }

  static std::size_t size() {
    return sizeof(block);
  }

  // Serialize/hash the full packed value (both 61-bit lanes), consistent with
  // size() == sizeof(block) and with the other x2 bundle types (fp107x2 /
  // fp2x128x2 / z2k64x2). Using sizeof(uint64_t) here would drop the high lane.
  template<typename IO>
  void send(IO *netio) {
    netio->send_data(&val, sizeof(block));
  }

  template<typename IO>
  void recv(IO *netio) {
    netio->recv_data(&val, sizeof(block));
  }

  block hash() {
    return Hash::hash_for_block(&val, sizeof(block));
  }


  // c = (a1 || a2) mod PR
  block vec_mod(block i) {
    i = _mm_add_epi64(
        (i & DoublePR()),
        _mm_srli_epi64(i, PR_bit_len));
    return vec_partial_mod(i);
  }

  // (c1 || c2) = (a1 || a2) (partial) mod PR
  block vec_partial_mod(const block i) const {
    return _mm_sub_epi64(
        i, _mm_andnot_si128(
        _mm_cmpgt_epi64(DoublePR(), i),
        DoublePR()));
  }
};

#endif // EMP_UTILS_FIELD_FP61x2_X2_H__
