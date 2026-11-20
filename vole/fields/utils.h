#ifndef VOLE_FIELDS_UTILS_H__
#define VOLE_FIELDS_UTILS_H__

#include "emp-tool/emp-tool.h"
#include <type_traits>
#include <vector>
#include <future>
#include <algorithm>
#include <cstring>

#if defined(__x86_64__) && defined(__BMI2__)
inline uint64_t vole_mul64(uint64_t a, uint64_t b, uint64_t *c) {
  return _mulx_u64((unsigned long long)a, (unsigned long long)b,
                   (unsigned long long *)c);
}
//
#else
inline uint64_t vole_mul64(uint64_t a, uint64_t b, uint64_t *c) {
  __uint128_t aa = a;
  __uint128_t bb = b;
  auto cc = aa * bb;
  *c = cc >> 64;
  return (uint64_t)cc;
}
#endif

// chi[i] = seed^(i+1) for i < 8, then chi[i] = chi[i-8] * seed^9, i.e.
//   chi[i] = seed^(i + 1 + floor(i/8)),
// distinct and increasing exponents up to ~1.125 * size (the multiplier is
// chi[7] * seed = seed^9, not seed^8; both parties use this function, so this
// is the coefficient vector of the protocol, and the polynomial-hash bound
// grows by the factor 1.125 only). Eight independent chains so the
// multiplications pipeline instead of forming one latency-bound dependency
// chain (emp-tool's block version does the same with four chains); ~4-5x
// faster for GF(2^128) gfmul.
template<typename T>
void field_uni_hash_coeff_gen(T *chi, T seed, int size) {
  if (size <= 0) return;
  chi[0] = seed;
  int lead = size < 8 ? size : 8;
  for (int i = 1; i < lead; ++i) chi[i] = chi[i - 1] * seed;
  if (size <= 8) return;
  T s8 = chi[7] * seed;  // seed^9 (see above)
  int i = 8;
  for (; i + 8 <= size; i += 8) {
    chi[i]     = chi[i - 8] * s8;
    chi[i + 1] = chi[i - 7] * s8;
    chi[i + 2] = chi[i - 6] * s8;
    chi[i + 3] = chi[i - 5] * s8;
    chi[i + 4] = chi[i - 4] * s8;
    chi[i + 5] = chi[i - 3] * s8;
    chi[i + 6] = chi[i - 2] * s8;
    chi[i + 7] = chi[i - 1] * s8;
  }
  for (; i < size; ++i) chi[i] = chi[i - 8] * s8;
}

// sum_{j < n} chi[j] * x(j) with chi = field_uni_hash_coeff_gen(seed, n),
// i.e. chi[j] = seed^(j + 1 + floor(j/8)), by Horner's rule instead (one
// multiplication per term instead of two, and no coefficient array). With
// j = 8m + r the exponent is 9m + (r + 1), so chain r collects
// x(8m + r) * (seed^9)^m and is weighted by seed^(r+1) at the end; eight
// interleaved chains keep the multiplier busy. Field arithmetic is exact, so
// the result equals the coefficient-array form bit for bit.
//
// Fields whose multiplication is long (FP192: ~80 instructions inlined) set
// poly_hash_horner<T> to false: eight such chains overflow the out-of-order
// window, and on Zen 5 the coefficient-array form measured ~5% faster there
// (the same value either way).
template<typename T> struct poly_hash_horner { static constexpr bool value = true; };

template<typename T, typename Get>
T field_poly_hash(T seed, std::size_t n, Get x) {
  if constexpr (!poly_hash_horner<T>::value) {
    std::vector<T> chi(n);
    field_uni_hash_coeff_gen(chi.data(), seed, (int)n);
    T res;
    res.setZero();
    for (std::size_t i = 0; i < n; ++i) res = res + chi[i] * x(i);
    return res;
  }
  T pw[8];  // pw[r] = seed^(r+1)
  pw[0] = seed;
  for (int r = 1; r < 8; ++r) pw[r] = pw[r - 1] * seed;
  const T s9 = pw[7] * seed;
  T acc[8];
  for (auto &a : acc) a.setZero();
  const std::size_t q = n / 8, rem = n % 8;
  for (std::size_t r = 0; r < rem; ++r) acc[r] = x(8 * q + r);
  for (std::size_t m = q; m-- > 0;)
    for (int r = 0; r < 8; ++r) acc[r] = acc[r] * s9 + x(8 * m + r);
  T res;
  res.setZero();
  for (int r = 0; r < 8; ++r) res = res + pw[r] * acc[r];
  return res;
}

// chi[e] = seed^(e + 1 + floor(e/8)) (field_uni_hash_coeff_gen) for
// 0 <= e < 2^bits, in a fixed sequence of operations: the exponent minus one,
// e + floor(e/8) < 2^(bits+1), is scanned over bits + 1 steps that each
// square and multiply; its bits only select (by index, not by branch) which
// product is kept. Used for the MPFSS receiver's coefficient at its secret
// point.
template<typename T>
T field_hash_coeff_fixed(T seed, std::size_t e, int bits) {
  const std::size_t x = e + (e >> 3);
  T r = seed, b = seed;
  for (int i = 0; i <= bits; ++i) {
    T opt[2] = {r, r * b};
    r = opt[(x >> i) & 1];
    b = b * b;
  }
  return r;
}

template<typename T>
T field_inn_prdt_sum_red(const T *a, const T *b, int size) {
  T res = a[0] * b[0];
  for(int i = 1; i < size; ++i) {
    res = res + (a[i] * b[i]);
  }
  return res;
}

// Restrict a freshly sampled Delta to the field/ring's admissible key space.
// No-op for the field types (any field element is a valid Delta); the Z_{2^k}
// ring overrides this (see z2k.h) to clamp Delta to its s-bit authentication
// sub-ring, because COPE only decomposes that many low bits of Delta.
template<typename T> inline void restrict_delta(T&) {}

// ---------------------------------------------------------------------------
// Ring (Z_{2^k}) vs. field dispatch.
//
// The malicious VOLE checks differ over a ring (MozZ2karella, Fig. 5) from those
// over a field (Wolverine): a field uses a single random-coefficient universal
// hash, but over Z_{2^k} the random ring coefficients have zero divisors, so the
// ring path instead runs (A) a GGM-tree tag check over GF(2^128) and (B) a
// linear-combination check with *binary* coefficients. `field_traits<T>::is_ring`
// gates these; it is specialized to `true` only for the ring type (see z2k.h).
// Call sites either dispatch on std::integral_constant<bool, is_ring>{} or use
// `if constexpr`, so the ring-only code is never instantiated for field types.
template<typename T> struct field_traits { static constexpr bool is_ring = false; };

// GGM leaf width (spfss_sender.h / spfss_recver.h). A leaf is one 128-bit block,
// mapped to a field element by from_block; that is a uniform element only when
// the field has at most ~128 bits. Fields wider than a block (FP192) set this to
// true and expand every leaf to two blocks (from_2blocks), otherwise the value
// at the punctured position, v[a] = w[a] - b*Delta, would be confined to
// [0, 2^128) and leak Delta. Specialized next to the field type.
template<typename T> struct wide_leaf { static constexpr bool value = false; };

// Primal-LPN matrix entries (lpn.h). Unit entries by default; a field set to
// true gives every entry its own pseudorandom coefficient of type
// lpn_coefficient<T>::type (built with from_block). GF(2^128) needs this: with
// 0/1 entries its LPN instance splits into 128 binary instances, one per bit
// plane, sharing one noise support. Prime fields have no such planes.
template<typename T> struct lpn_coefficients { static constexpr bool value = false; };
template<typename T> struct lpn_coefficient { using type = T; };
template<typename T> using lpn_coefficient_t = typename lpn_coefficient<T>::type;

// Uniform bits per commitment-matrix coefficient, as used by the Theorem-1
// binding bound (cvole.h: cvole_n_com_for). Defaults to log|F|; a field whose
// coefficients carry less entropy than its size (FP192: one PRP block) says so
// next to its type.
template<typename T> struct com_coef_bits { static constexpr double value = (double)T::PR_bit_len; };

// Coefficients of the base-sVOLE universal-hash check (base_svole.h). Over a
// field: the polynomial hash chi_j = seed^(j+1). Over Z_{2^k} that hash
// degenerates: an even seed has seed^128 = 0, so every coefficient from j = 127
// on vanishes and those pairs go unchecked (with probability 1/2 per run). The
// ring therefore draws independent uniform ring elements from PRG(seed): for
// the first nonzero error e_j with v trailing zeros the combination is uniform
// on 2^v Z_{2^k}, so an error in the low k - v bits passes with probability
// 2^-(k-v) (the SPDZ2k bound; the high s bits are not authenticated).
template<typename T>
void check_coeff_gen(T *chi, T seed, int size) {
  if constexpr (field_traits<T>::is_ring) {
    static_assert(sizeof(T) >= sizeof(emp::block), "ring seed must fill a block");
    emp::block b;
    std::memcpy(&b, &seed, sizeof(b));
    emp::PRG prg(&b);
    prg.random_data_unaligned(chi, size * (int)sizeof(T));
  } else {
    field_uni_hash_coeff_gen(chi, seed, size);
  }
}

// Protocol hooks keyed on the value||MAC bundle type FPS. Defaults describe a
// full-field value: every MPFSS point consumes one base pair (its value), the
// consistency-check mask is one base pair, base values are uniform field
// elements, and base keys need no normalisation. The F_2 bundle (fields/f2.h)
// specialises all four to obtain Ferret's correlated OT.
template<typename FPS> struct vole_traits {
  static constexpr bool unit_point_value = false;   // point value = base pair's value
  static constexpr std::size_t mask_pairs = 1;      // base pairs forming the check mask
  using mask_t = FPS;
  static mask_t pack_mask(const FPS *p) { return p[0]; }
  template <typename FP> static FP pack_mask_key(const FP *k) { return k[0]; }
  template <typename FP, typename PRNG> static void rand_base_value(FP &x, PRNG &prg) { x.rand(prg); }
  template <typename FP> static void normalize_base_keys(FP *, std::size_t) {}
};

// Fixed key for the GGM leaf-tag PRG used by the ring GGM-tree check. The tag map
// is a fixed-key AES permutation (injective), giving the right-half-injectivity
// the check requires (MozZ2karella Def. 3 / Thm. 5).
const static emp::block GGM_TAG_KEY =
    emp::makeBlock(0x9e3779b97f4a7c15ULL, 0xbf58476d1ce4e5b9ULL);

// Tag a batch of GGM leaves in place: t_j <- AES_{GGM_TAG_KEY}(leaf_j).
inline void ggm_leaf_tag(emp::block *tag, const emp::block *leaf, int n) {
  emp::PRP prp(GGM_TAG_KEY);
  for (int i = 0; i < n; ++i) tag[i] = leaf[i];
  prp.permute_block(tag, n);
}

// Derive a distinct per-tree seed from a shared per-thread seed: AES_seed(idx).
// Both parties hold the same `seed` and tree index, so they derive the same
// binary chi (ocelot derives per-instance seeds the same way, prover.rs:247).
inline emp::block ring_tree_seed(emp::block seed, std::size_t idx) {
  emp::PRP prp(seed);
  emp::block t = emp::makeBlock(0, (uint64_t)idx);
  prp.permute_block(&t, 1);
  return t;
}

// Expand `seed` to a {0,1}^n vector of Hamming weight n/2 by sampling n/2 distinct
// indices (port of ocelot spvole/prover.rs:190-206). Deterministic in `seed`, so
// both parties derive the identical vector from the shared per-tree seed.
inline void gen_binary_chi(bool *chi, emp::block seed, int n) {
  for (int i = 0; i < n; ++i) chi[i] = false;
  emp::PRG prg(&seed);
  int set = 0;
  while (set < n / 2) {
    uint32_t r;
    prg.random_data_unaligned(&r, sizeof(r));
    uint32_t idx = r % (uint32_t)n;
    if (!chi[idx]) { chi[idx] = true; ++set; }
  }
}

#endif // VOLE_FIELDS_UTILS_H__
