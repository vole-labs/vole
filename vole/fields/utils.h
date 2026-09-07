#ifndef VOLE_FIELDS_UTILS_H__
#define VOLE_FIELDS_UTILS_H__

#include "emp-tool/emp-tool.h"
#include <type_traits>
#include <vector>
#include <future>
#include <algorithm>
#include <cstring>

#if defined(__x86_64__) && defined(__BMI2__)
inline uint64_t mul64(uint64_t a, uint64_t b, uint64_t *c) {
  return _mulx_u64((unsigned long long)a, (unsigned long long)b,
                   (unsigned long long *)c);
}
//
#else
inline uint64_t mul64(uint64_t a, uint64_t b, uint64_t *c) {
  __uint128_t aa = a;
  __uint128_t bb = b;
  auto cc = aa * bb;
  *c = cc >> 64;
  return (uint64_t)cc;
}
#endif

// chi[i] = seed^(i+1). Computed as eight independent chains
// (chi[i+8] = chi[i] * seed^8) so the multiplications pipeline instead of
// forming one latency-bound dependency chain (emp-tool's block version does
// the same with four chains); ~4-5x faster for GF(2^128) gfmul.
template<typename T>
void field_uni_hash_coeff_gen(T *chi, T seed, int size) {
  if (size <= 0) return;
  chi[0] = seed;
  int lead = size < 8 ? size : 8;
  for (int i = 1; i < lead; ++i) chi[i] = chi[i - 1] * seed;
  if (size <= 8) return;
  T s8 = chi[7] * seed;  // seed^8
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
// We build for C++11 (no `if constexpr`), so the call sites dispatch on
// std::integral_constant<bool, field_traits<FP>::is_ring>{} and the ring-only
// code is never instantiated for the field types.
template<typename T> struct field_traits { static constexpr bool is_ring = false; };

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
