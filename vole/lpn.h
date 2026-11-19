#ifndef _LPN_FP_H__
#define _LPN_FP_H__

#include "emp-tool/emp-tool.h"
#include "vole/fields/utils.h"
#include <string>
#include <string_view>

// Primal-LPN expansion: out[o] += sum_{j < d} a_{o,j} * pre[idx_{o,j}] with
// d = 10 pseudorandom indices into the k-entry secret `pre`.
//
// The matrix is public and comes from a 128-bit root agreed per VOLE instance
// (VoleTriple::set_lpn_seed / its setup exchange). Every compute() call is one
// round and reads fresh columns: round r of a stage uses the AES key
// RO(domain, root || r), so stages (distinct domains), rounds, and sessions
// (distinct roots) never share a matrix.
//
// k must be a power of two, so idx = r & (k - 1) is uniform on [0, k).
//
// Entries a_{o,j} are 1, except over GF(2^128) (lpn_coefficients<FP>), where
// each entry is its own pseudorandom field element: with 0/1 entries the
// instance would split into 128 binary instances, one per bit plane, sharing
// one noise support. The map is linear, so MAC = key + value * Delta carries
// over when both lanes are multiplied by the same a_{o,j}.
template <int d = 10>
class LpnFp {
public:
  int k, n;
  ThreadPool *pool;
  int threads;
  block root;
  std::string domain;
  uint64_t round_ = 0;  // compute() calls so far
  uint32_t k_mask;

  LpnFp(int n, int k, ThreadPool *pool, int threads, block root,
        std::string_view domain)
      : k(k), n(n), pool(pool), threads(threads), root(root), domain(domain) {
    if (k <= 0 || (k & (k - 1)) != 0)
      error("LpnFp: k must be a power of two");
    k_mask = (uint32_t)k - 1;
  }

  // Outputs are processed in batches of B = 32; batch b draws its 320 indices
  // from one PRP call on the 80 blocks makeBlock(b, m), m < 80, and (over
  // GF(2^128)) its 320 coefficients from a second PRP on makeBlock(b, m),
  // m < 320. Indices depend only on (round, b, position in batch), so the
  // matrix is independent of the thread count; threads partition batches.
  static constexpr int B = 32;
  static constexpr int NB = (B * d + 3) / 4;  // 80 index blocks per batch
  static constexpr int NC = B * d;            // 320 coefficient blocks per batch

  // AES keys of round `r`: matrix indices and (GF(2^128)) coefficients.
  block round_key(uint64_t r) const {
    return RO(domain, root).absorb(std::string_view("index")).absorb(r).squeeze_block();
  }
  block coef_key(uint64_t r) const {
    return RO(domain, root).absorb(std::string_view("coefficient")).absorb(r).squeeze_block();
  }

  template <typename FP>
  void compute_batch(FP *K, const FP *preK, int64_t b, int cnt, PRP *prp,
                     PRP *cprp, block *tmp, block *coef) {
    for (int m = 0; m < NB; ++m) tmp[m] = makeBlock((uint64_t)b, (uint64_t)m);
    prp->permute_block(tmp, NB);
    const uint32_t *r = (const uint32_t *)(tmp);
    const int64_t base = b * (int64_t)B;
    if constexpr (lpn_coefficients<FP>::value) {
      for (int m = 0; m < NC; ++m) coef[m] = makeBlock((uint64_t)b, (uint64_t)m);
      cprp->permute_block(coef, NC);
      for (int o = 0; o < cnt; ++o) {
        FP acc = K[base + o];
        for (int j = 0; j < d; ++j) {
          lpn_coefficient_t<FP> a;
          a.from_block(coef[o * d + j]);
          acc.add_raw(preK[r[o * d + j] & k_mask] * a);
        }
        K[base + o] = acc;
      }
    } else if constexpr (requires { typename FP::Acc; }) {
      (void)cprp; (void)coef;
      // Element types without spare bits (FP192) supply a wider lazy
      // accumulator: d unreduced adds, one reduction per output.
      for (int o = 0; o < cnt; ++o) {
        typename FP::Acc acc(K[base + o]);
#pragma GCC unroll 10
        for (int j = 0; j < d; ++j) acc.add(preK[r[o * d + j] & k_mask]);
        K[base + o] = acc.get();
      }
    } else {
      (void)cprp; (void)coef;
      for (int o = 0; o < cnt; ++o) {
        FP acc = K[base + o];
#pragma GCC unroll 10
        for (int j = 0; j < d; ++j) {
          acc.add_raw(preK[r[o * d + j] & k_mask]);
          if ((j + 1) % (int)FP::lazy_adds == 0) acc.reduce();
        }
        if (d % (int)FP::lazy_adds != 0) acc.reduce();
        K[base + o] = acc;
      }
    }
  }

  // Batches [b0, b1) of this thread, round keys given.
  template <typename FP>
  void task(FP *K, const FP *preK, int64_t b0, int64_t b1, block ik, block ck) {
    PRP prp(ik);
    PRP cprp(ck);
    alignas(16) block tmp[NB];
    alignas(16) block coef[lpn_coefficients<FP>::value ? NC : 1];
    for (int64_t b = b0; b < b1; ++b) {
      int cnt = (int)std::min<int64_t>(B, (int64_t)n - b * B);
      compute_batch<FP>(K, preK, b, cnt, &prp, &cprp, tmp, coef);
    }
  }

  // One round: K[0, n) += A_round * preK[0, k).
  template <typename FP>
  void compute(FP *K, const FP *preK) {
    const block ik = round_key(round_);
    const block ck = lpn_coefficients<FP>::value ? coef_key(round_) : zero_block;
    ++round_;
    std::vector<std::future<void>> fut;
    int64_t nbatch = ((int64_t)n + B - 1) / B;
    int64_t width = (nbatch + threads - 1) / threads;
    for (int i = 0; i < threads - 1; ++i) {
      int64_t b0 = i * width, b1 = std::min<int64_t>((i + 1) * width, nbatch);
      if (b0 >= b1) break;
      fut.push_back(pool->enqueue([this, K, preK, b0, b1, ik, ck]() { task(K, preK, b0, b1, ik, ck); }));
    }
    int64_t b0 = (int64_t)(threads - 1) * width;
    if (b0 < nbatch) task(K, preK, b0, nbatch, ik, ck);
    for (auto &f : fut)
      f.get();
  }
};
#endif
