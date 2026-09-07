#ifndef _LPN_FP_H__
#define _LPN_FP_H__

#include "emp-tool/emp-tool.h"

template <int d = 10>
class LpnFp {
public:
  int k, n;
  ThreadPool *pool;
  int threads;
  block seed;

  uint32_t k_mask;
  LpnFp(int n, int k, ThreadPool *pool, int threads, block seed = zero_block) {
    this->k = k;
    this->n = n;
    this->pool = pool;
    this->threads = threads;
    this->seed = seed;

    // Smallest 2^m - 1 >= k; for a power-of-two k exactly k - 1, so that
    // r & k_mask is already uniform on [0, k) and no fold-back is needed.
    k_mask = 1;
    while (k_mask < (uint32_t)k) {
      k_mask <<= 1;
      k_mask = k_mask | 0x1;
    }
    if ((k & (k - 1)) == 0) k_mask = (uint32_t)k - 1;
  }

  // Index generation: outputs are processed in batches of B = 32; batch b
  // draws its 320 indices from one PRP call on the 80 blocks
  // makeBlock(b, m), m < 80 (a batch large enough for emp-tool's wide AES
  // tiles). Indices depend only on (b, position in batch), so the LPN matrix
  // is independent of the thread count on either party; threads partition
  // batches, never split one. Each index is r & k_mask, then r >= k ? r - k
  // (skipped when k is a power of two).
  static constexpr int B = 32;
  static constexpr int NB = (B * 10 + 3) / 4;  // 80 blocks per batch

  template<typename FP, bool POW2>
  void compute_batch(FP *K, const FP *preK, int64_t b, int cnt, PRP *prp, block *tmp) {
    for (int m = 0; m < NB; ++m) tmp[m] = makeBlock((uint64_t)b, (uint64_t)m);
    prp->permute_block(tmp, NB);
    const uint32_t *r = (const uint32_t *)(tmp);
    const int64_t base = b * (int64_t)B;
    for (int o = 0; o < cnt; ++o) {
      FP acc = K[base + o];
#pragma GCC unroll 10
      for (int j = 0; j < 10; ++j) {
        int idx = (int)(r[o * 10 + j] & k_mask);
        if (!POW2) idx = (idx >= k) ? idx - k : idx;
        acc.add_raw(preK[idx]);
        if ((j + 1) % (int)FP::lazy_adds == 0) acc.reduce();
      }
      if (10 % (int)FP::lazy_adds != 0) acc.reduce();
      K[base + o] = acc;
    }
  }

  // Batches [b0, b1) of this thread.
  template<typename FP>
  void task(FP *K, const FP *preK, int64_t b0, int64_t b1) {
    PRP prp(seed);
    alignas(16) block tmp[NB];
    const bool pow2 = (k & (k - 1)) == 0;
    for (int64_t b = b0; b < b1; ++b) {
      int cnt = (int)std::min<int64_t>(B, (int64_t)n - b * B);
      if (pow2) compute_batch<FP, true>(K, preK, b, cnt, &prp, tmp);
      else      compute_batch<FP, false>(K, preK, b, cnt, &prp, tmp);
    }
  }

  template<typename FP>
  void compute(FP *K, const FP* preK) {
    std::vector<std::future<void>> fut;
    int64_t nbatch = ((int64_t)n + B - 1) / B;
    int64_t width = (nbatch + threads - 1) / threads;
    for (int i = 0; i < threads - 1; ++i) {
      int64_t b0 = i * width, b1 = std::min<int64_t>((i + 1) * width, nbatch);
      if (b0 >= b1) break;
      fut.push_back(pool->enqueue([this, K, preK, b0, b1]() { task(K, preK, b0, b1); }));
    }
    int64_t b0 = (int64_t)(threads - 1) * width;
    if (b0 < nbatch) task(K, preK, b0, nbatch);
    for (auto &f : fut)
      f.get();
  }

};
#endif
