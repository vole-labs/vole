#ifndef EMP_UTILS_FIELD_FP61_H__
#define EMP_UTILS_FIELD_FP61_H__

#include "emp-tool/emp-tool.h"
#include "vole/fields/utils.h"
#include <cstdint>

using namespace emp;

// finite field Fp, p = 2^61 - 1
class FP61 {
public:
  static constexpr uint64_t PR_bit_len       = 61;
  static constexpr uint64_t slot_stride_bits = 61;   // = PR_bit_len
  static constexpr uint64_t PR               = 2305843009213693951;  // 2^61 - 1
  static constexpr uint64_t PR_mask          = (1ULL << 61) - 1;
  static constexpr uint64_t PR_num_pack      = 1ULL;
  static constexpr uint64_t DS_byte_len      = 8;
  // A `block` is not a literal type, so DoublePR can't be `constexpr`. Expose it
  // through an inline accessor (implicitly inline -> single definition across
  // translation units) instead of a file-scope static, which would otherwise
  // produce duplicate symbols when two .cpp files include this header.
  static block DoublePR() {
    return makeBlock(2305843009213693951, 2305843009213693951);
  }

  uint64_t val;

  FP61() { val = 0; }

  FP61(uint64_t input, bool mod_it = true) {
    if (mod_it)
      val = mod(input);
    else
      val = input;
  }

  void operator=(const uint64_t rhs) { this->val = mod(rhs); }

  void operator=(const FP61 rhs) { this->val = rhs.val; }

  void setZero() { this->val = 0ULL; }

  void from_block(block rhs) {
    uint64_t x = _mm_extract_epi64(rhs, 0);
    this->val = mod(x);
  }

  void assign_no_mod(const uint64_t rhs) { this->val = rhs; }

  uint64_t value() {
    return val;
  }

  void setHigh(const uint64_t rhs) {
    this->val = rhs;
  }

  void setHigh(const FP61 &rhs) {
    this->val = rhs.val;
  }

  FP61 getHigh() {
    return FP61(this->val, false);
  }

  void setLow(const FP61 &rhs) {
    this->val = rhs.val;
  }

  FP61 getLow() {
    return FP61(this->val, false);
  }

  bool operator==(const uint64_t rhs) const {
    return (this->val == mod(rhs)) ? true : false;
  }

  bool operator==(const FP61 rhs) const {
    return (this->val == rhs.val) ? true : false;
  }

  bool operator!=(const uint64_t rhs) const {
    return (this->val == mod(rhs)) ? false : true;
  }

  bool operator!=(const FP61 rhs) const {
    return (this->val == rhs.val) ? false : true;
  }

  FP61 operator+(const FP61 rhs) const {
    FP61 res(*this);
    res.val = add_mod(res.val, rhs.val);
    return res;
  }

  FP61 operator+(const uint64_t rhs) const {
    FP61 res(*this);
    res.val = add_mod(res.val, rhs);
    return res;
  }

  FP61 operator-(const FP61 rhs) const {
    FP61 res(*this);
    res.val = add_mod(res.val, PR - rhs.val);
    return res;
  }

  FP61 operator-(const uint64_t rhs) const {
    FP61 res(*this);
    res.val = add_mod(res.val, PR - rhs);
    return res;
  }

  FP61 operator*(const FP61 rhs) const {
    FP61 res(*this);
    res.val = mult_mod(res.val, rhs.val);
    return res;
  }

  FP61 operator*(const uint64_t rhs) const {
    FP61 res(*this);
    res.val = mult_mod(res.val, rhs);
    return res;
  }

  FP61 operator<<(const int n_pos) const {
    FP61 res(*this);
    res.val = mult_mod(res.val, (1ULL<<n_pos));
    return res;
  }

  // Lazy reduction: up to `lazy_adds` add_raw() calls on a reduced value stay
  // below 2^64 (6 * 2^61 < 2^64); reduce() brings the sum back below p.
  static constexpr unsigned lazy_adds = 5;
  void add_raw(const FP61 &rhs) { val += rhs.val; }
  void reduce() { val = mod(val); }

  FP61 negate() const {
    FP61 res(*this);
    if(res.val != 0)
      res.val = PR - res.val;
    return res;
  }

  template<typename IO>
  void send(IO *netio) {
    netio->send_data(&val, sizeof(uint64_t));
  }

  template<typename IO>
  void recv(IO *netio) {
    netio->recv_data(&val, sizeof(uint64_t));
  }

  template<typename PRNG>
  void rand(PRNG &prg) {
    uint64_t raw;
    prg.random_data_unaligned(&raw, sizeof(uint64_t));
    raw = mod(raw);
    this->val = raw;
  }

  block hash() {
    return Hash::hash_for_block(&val, sizeof(uint64_t));
  }

  FP61 inv() const {
    FP61 res(*this);
    res.val = mod_inv(res.val);
    return res;
  }

  uint64_t mod(const uint64_t x) const {
    uint64_t i = (x & PR) + (x >> PR_bit_len);
    return (i >= PR) ? i - PR : i;
  }

  // c = val + a mod PR
  uint64_t add_mod(const uint64_t a, const uint64_t b) const {
    uint64_t res = a + b;
    return (res >= PR) ? (res - PR) : res;
  }

  // c = a * b mod PR
  uint64_t mult_mod(const uint64_t a, const uint64_t b) const {
    uint64_t c = 0;
    uint64_t e = vole_mul64(a, b, (uint64_t *)&c);
    uint64_t res = (e & PR) + ((e >> PR_bit_len) ^ (c << (64 - PR_bit_len)));
    return (res >= PR) ? (res - PR) : res;
  }

  // from SCL
  uint64_t mod_inv(const uint64_t a) const {
    if(a == 0)
      error("zero know invertible");

    int64_t k = 0;
    int64_t new_k = 1;
    int64_t r = (int64_t)PR;
    int64_t new_r = (int64_t)a;

    while(new_r != 0) {
      int64_t q = r / new_r;
      // assign(k, new_k, q)
      int64_t tmp = new_k;
      new_k = k - q * tmp;
      k = tmp;
      // assign(r, new_r, q)
      tmp = new_r;
      new_r = r - q * tmp;
      r = tmp;
    }

    if(k < 0) k = k + (int64_t)PR;
    return (uint64_t)k;
  }

  static std::size_t size() {
    return sizeof(uint64_t);
  }

  static uint64_t copy_compose(uint64_t basis) {
    return basis;
  }


  // c = (a1 * b || a2 * b) mod PR
  static block mult_mod(block a, uint64_t b) {
    uint64_t H = _mm_extract_epi64(a, 1);
    uint64_t L = _mm_extract_epi64(a, 0);
    block bs[2];
    uint64_t *is = (uint64_t *)(bs);
    is[1] = vole_mul64(H, b, (uint64_t *)(is + 3));
    is[0] = vole_mul64(L, b, (uint64_t *)(is + 2));
    block t1 = bs[0] & DoublePR();
    block t2 = _mm_srli_epi64(bs[0], PR_bit_len) ^
               _mm_slli_epi64(bs[1], 64 - PR_bit_len);
    block res = _mm_add_epi64(t1, t2);
    return vec_partial_mod(res);
  }

  // c = (a1 + b || a2 + b) mod PR
  static block add_mod(block a, uint64_t b) {
    block res = _mm_add_epi64(a, _mm_set_epi64((__m64)b, (__m64)b));
    return vec_partial_mod(res);
  }

  // c = (a1 + b1 || a2 + b2) mod PR
  static block add_mod(block a, block b) {
    block res = _mm_add_epi64(a, b);
    return vec_partial_mod(res);
  }

  // c = (a1 || a2) mod PR
  static block vec_mod(block i) {
    i = _mm_add_epi64((i & DoublePR()), _mm_srli_epi64(i, PR_bit_len));
    return vec_partial_mod(i);
  }

  // (c1 || c2) = (a1 || a2) (partial) mod PR
  static block vec_partial_mod(block i) {
    return _mm_sub_epi64(
        i, _mm_andnot_si128(_mm_cmpgt_epi64(DoublePR(), i), DoublePR()));
  }
};

#endif // EMP_UTILS_FIELD_FP61_H__
