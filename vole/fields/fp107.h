#ifndef EMP_UTILS_FIELD_FP107_H__
#define EMP_UTILS_FIELD_FP107_H__

// Finite field F_p, p = 2^107 - 1 (the 14th Mersenne prime, M107).
//
// Arithmetic backbone:
//   - Storage: one `unsigned __int128` per field element (107 useful bits;
//              top 21 bits are always zero for a reduced value).
//   - Add/sub: `__int128` arithmetic + a single conditional subtract of p.
//   - Mul:     a * b is up to 214 bits. We compute the four 64x64 -> 128
//              partial products, assemble them into a 256-bit (high, low)
//              pair (held as two `__int128`s), then split into two
//              107-bit chunks c0, c1 and reduce c0 + c1 mod p with at
//              most two conditional subtracts.
//   - Inverse: extended Euclidean on `__int128` with sign tracking.
//
// Storage / wire layout:
//   `val` holds the 107-bit value packed into a `__uint128_t`, with
//   bits 107..127 always zero. Helpers `lo64()` / `hi64()` extract the
//   low and (43-bit) high halves for serialization.

#include "emp-tool/emp-tool.h"
#include "vole/fields/utils.h"
#include <cstdint>
#include <cstring>

using namespace emp;

class FP107 {
 public:
    using u128 = unsigned __int128;
    using i128 = __int128;

    // 2^107 - 1.
    static constexpr uint64_t PR_bit_len       = 107;
    static constexpr uint64_t slot_stride_bits = 107;            // = PR_bit_len
    static constexpr u128     PR               = (((u128)1) << 107) - 1;
    static constexpr u128     HALF_PR    = ((((u128)1) << 107) - 1) / 2;
    static constexpr u128     PR_mask    = (((u128)1) << 107) - 1;
    static constexpr uint64_t PR_num_pack = 1;
    static constexpr uint64_t DS_byte_len = 16;

    u128 val;

    FP107() : val(0) {}
    FP107(u128 v, bool reduce = true) : val(reduce ? mod(v) : v) {}
    // Non-explicit so int / uint64_t literals convert implicitly (matching
    // FP61's surface). Always reduces.
    FP107(uint64_t v)                : val((u128)v) {}
    FP107(int v)                     : val(v < 0 ? PR - ((u128)(-v) % PR) : (u128)v) {}
    // Explicit (lo, hi) split-decoder factory — explicit so it doesn't
    // collide with FP107(uint64_t, bool=true) overload resolution.
    static FP107 from_lo_hi(uint64_t lo, uint64_t hi) {
        FP107 r; r.val = mod(((u128)hi << 64) | (u128)lo); return r;
    }

    void operator=(const FP107& rhs) { val = rhs.val; }
    void operator=(uint64_t rhs)      { val = mod((u128)rhs); }

    void setZero() { val = 0; }

    void assign_no_mod(const u128 rhs) { val = rhs; }

    u128 value() const { return val; }

    // SPDZ MAC layout helpers: a single FP107 holds either the value or the
    // MAC -- the (value, MAC) pair lives in FP107x2. setHigh/getHigh exist
    // for symmetry with FP61 so the same templates compile.
    void setHigh(const u128 rhs) { val = rhs; }
    void setHigh(const FP107 &rhs) { val = rhs.val; }

    FP107 getHigh() {
        FP107 r; r.val = val; return r;
    }

    void setLow(const FP107 &rhs) { val = rhs.val; }

    FP107 getLow() {
        FP107 r; r.val = val; return r;
    }

    bool operator==(const FP107& rhs) const { return val == rhs.val; }
    bool operator!=(const FP107& rhs) const { return val != rhs.val; }
    bool operator==(const u128& rhs) const { return val == mod(rhs); }
    bool operator!=(const u128& rhs) const { return val != mod(rhs); }

    FP107 operator+(const FP107& rhs) const {
        u128 s = val + rhs.val;
        if (s >= PR) s -= PR;
        FP107 r; r.val = s; return r;
    }
    FP107 operator-(const FP107& rhs) const {
        FP107 r;
        r.val = (val >= rhs.val) ? (val - rhs.val) : (val + PR - rhs.val);
        return r;
    }
    FP107 operator-() const {
        FP107 r;
        r.val = (val == 0) ? 0 : (PR - val);
        return r;
    }
    FP107 operator*(const FP107& rhs) const {
        FP107 r; r.val = mult_mod(val, rhs.val); return r;
    }

    // Scalar overloads (uint64_t / u128) — let cope/spfss code that builds
    // values from raw uint64 bit operations work without explicit casts.
    FP107 operator+(uint64_t rhs) const { return *this + FP107(rhs); }
    FP107 operator-(uint64_t rhs) const { return *this - FP107(rhs); }
    FP107 operator*(uint64_t rhs) const { FP107 r; r.val = mult_mod(val, (u128)rhs); return r; }
    FP107 operator+(u128 rhs) const { FP107 r; r.val = add_mod(val, mod(rhs)); return r; }
    FP107 operator-(u128 rhs) const { FP107 r; r.val = sub_mod(val, mod(rhs)); return r; }
    FP107 operator*(u128 rhs) const { FP107 r; r.val = mult_mod(val, mod(rhs)); return r; }

    FP107 operator<<(const int n_pos) const {
        FP107 r; r.val = mult_mod(val, ((u128)1) << n_pos); return r;
    }

    // Comparison with raw scalars (matches FP61).
    bool operator==(uint64_t rhs) const { return val == mod((u128)rhs); }
    bool operator!=(uint64_t rhs) const { return val != mod((u128)rhs); }

    FP107& operator+=(const FP107& rhs) { *this = *this + rhs; return *this; }
    FP107& operator-=(const FP107& rhs) { *this = *this - rhs; return *this; }
    FP107& operator*=(const FP107& rhs) { *this = *this * rhs; return *this; }

    // Lazy reduction: values are < 2^107 in a 128-bit word, so 2^21 - 2 raw
    // adds fit; reduce() is the single Mersenne fold.
    static constexpr unsigned lazy_adds = 1024;
    void add_raw(const FP107 &rhs) { val += rhs.val; }
    void reduce() { val = mod(val); }

    FP107 negate() const {
        FP107 r;
        r.val = (val == 0) ? 0 : (PR - val);
        return r;
    }

    // Multiplicative inverse using extended Euclidean. Returns 0 for val == 0.
    FP107 inv() const {
        FP107 r; r.val = mod_inv(val); return r;
    }

    template<typename PRNG>
    void rand(PRNG& prg) {
        // Pull 16 bytes, mask to 107 bits, conditionally reduce.
        uint64_t lo, hi;
        prg.random_data(&lo, sizeof(uint64_t));
        prg.random_data(&hi, sizeof(uint64_t));
        u128 x = ((u128)hi << 64) | (u128)lo;
        x &= PR;             // bring into [0, 2^107 - 1]
        if (x == PR) x = 0;  // collapse the single non-canonical value
        val = x;
    }

    template<typename IO>
    void send(IO *netio) {
        netio->send_data(&val, sizeof(u128));
    }

    template<typename IO>
    void recv(IO *netio) {
        netio->recv_data(&val, sizeof(u128));
    }

    block hash() {
        return Hash::hash_for_block(&val, sizeof(u128));
    }

    void from_block(block rhs) {
        uint64_t lo = _mm_extract_epi64(rhs, 0);
        uint64_t hi = _mm_extract_epi64(rhs, 1);
        u128 x = ((u128)hi << 64) | (u128)lo;
        val = mod(x);
    }

    static u128 copy_compose(u128 basis) { return basis; }
    static uint64_t copy_compose(uint64_t basis) { return basis; }

    // Helpers: (lo, hi) decomposition for serialization. lo is 64 bits;
    // hi is the upper 43 bits, padded into a uint64_t.
    uint64_t lo64() const { return (uint64_t)(val & 0xFFFFFFFFFFFFFFFFULL); }
    uint64_t hi64() const { return (uint64_t)(val >> 64); }

    // Mersenne reduction: x mod (2^107 - 1) for x < 2^128 (a single u128).
    static u128 mod(u128 x) {
        u128 r = (x & PR) + (x >> 107);
        if (r >= PR) r -= PR;
        return r;
    }

    static u128 add_mod(u128 a, u128 b) {
        u128 s = a + b;
        return (s >= PR) ? (s - PR) : s;
    }
    static u128 sub_mod(u128 a, u128 b) {
        return (a >= b) ? (a - b) : (a + PR - b);
    }

    // Multiply two values in [0, 2^107) and reduce mod p.
    static u128 mult_mod(u128 a, u128 b) {
        uint64_t a_lo = (uint64_t)a;
        uint64_t a_hi = (uint64_t)(a >> 64);
        uint64_t b_lo = (uint64_t)b;
        uint64_t b_hi = (uint64_t)(b >> 64);

        u128 ll = (u128)a_lo * (u128)b_lo;       // bits   0..127
        u128 lh = (u128)a_lo * (u128)b_hi;       // bits  64..170
        u128 hl = (u128)a_hi * (u128)b_lo;       // bits  64..170
        u128 hh = (u128)a_hi * (u128)b_hi;       // bits 128..213

        u128 mid     = lh + hl;
        u128 mid_carry = (mid < lh) ? ((u128)1 << 64) : 0;

        u128 low_add = mid << 64;
        u128 low     = ll + low_add;
        u128 carry   = (low < ll) ? 1 : 0;

        u128 high = hh + (mid >> 64) + mid_carry + carry;

        u128 c0 = low & PR;
        u128 c1 = (low >> 107) | (high << 21);

        u128 s = c0 + c1;
        if (s >= PR) s -= PR;
        if (s >= PR) s -= PR;
        return s;
    }

    // Modular inverse via extended Euclidean over signed 128-bit ints.
    static u128 mod_inv(u128 a) {
        if (a == 0) return 0;

        i128 r_old = (i128)PR;
        i128 r_new = (i128)a;
        i128 k_old = 0;
        i128 k_new = 1;

        while (r_new != 0) {
            i128 q = r_old / r_new;
            i128 t;
            t = r_old - q * r_new; r_old = r_new; r_new = t;
            t = k_old - q * k_new; k_old = k_new; k_new = t;
        }

        if (k_old < 0) k_old += (i128)PR;
        return (u128)k_old;
    }


    static std::size_t size() { return sizeof(u128); }
};

#endif  // EMP_UTILS_FIELD_FP107_H__
