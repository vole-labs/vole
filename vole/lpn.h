#ifndef _LPN_FP_H__
#define _LPN_FP_H__

#include "emp-tool/emp-tool.h"
#include "vole/aes_ni.h"

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

    k_mask = 1;
    while (k_mask < (uint32_t)k) {
      k_mask <<= 1;
      k_mask = k_mask | 0x1;
    }
  }

  // Map a raw 32-bit PRP word to an index in [0, k).
  inline int fix_index(uint32_t x) const {
    int v = (int)(x & k_mask);
    return (v >= k) ? v - k : v;
  }

  template<typename FP>
  void __compute4(FP *K, const FP* preK, int i, PRP *prp) {
    block tmp[10];
    for (int m = 0; m < 10; ++m)
      tmp[m] = makeBlock(i, m);
    vole::aes_ecb_encrypt_blks(tmp, 10, &prp->aes);
    // Indices are masked as they are consumed. The reduction schedule folds
    // at compile time: raw adds, one reduce every FP::lazy_adds terms (2 per
    // output for fp61, matching the hand-unrolled loop in emp-zk).
    const uint32_t *r = (const uint32_t *)(tmp);
    FP tmpv[4];
    tmpv[0] = K[i];
    tmpv[1] = K[i+1];
    tmpv[2] = K[i+2];
    tmpv[3] = K[i+3];
#pragma GCC unroll 10
    for (int j = 0; j < 10; ++j) {
      tmpv[0].add_raw(preK[fix_index(r[4 * j + 0])]);
      tmpv[1].add_raw(preK[fix_index(r[4 * j + 1])]);
      tmpv[2].add_raw(preK[fix_index(r[4 * j + 2])]);
      tmpv[3].add_raw(preK[fix_index(r[4 * j + 3])]);
      if ((j + 1) % (int)FP::lazy_adds == 0) {
        tmpv[0].reduce(); tmpv[1].reduce(); tmpv[2].reduce(); tmpv[3].reduce();
      }
    }
    if (10 % (int)FP::lazy_adds != 0) {
      tmpv[0].reduce(); tmpv[1].reduce(); tmpv[2].reduce(); tmpv[3].reduce();
    }
    K[i] = tmpv[0];
    K[i+1] = tmpv[1];
    K[i+2] = tmpv[2];
    K[i+3] = tmpv[3];
  }

  template<typename FP>
  void __compute1(FP *K, const FP* preK, int i, PRP *prp) {
    block tmp[3];
    for (int m = 0; m < 3; ++m)
      tmp[m] = makeBlock(i, m);
    vole::aes_ecb_encrypt_blks(tmp, 3, &prp->aes);
    const uint32_t *r = (const uint32_t *)(tmp);
    FP tmpv = K[i];
#pragma GCC unroll 10
    for (int j = 0; j < 10; ++j) {
      tmpv.add_raw(preK[fix_index(r[j])]);
      if ((j + 1) % (int)FP::lazy_adds == 0) tmpv.reduce();
    }
    if (10 % (int)FP::lazy_adds != 0) tmpv.reduce();
    K[i] = tmpv;
  }

  template<typename FP>
  void task(FP *K, const FP* preK, int start, int end) {
    PRP prp(seed);
    int j = start;
    for (; j < end - 4; j += 4)
      __compute4(K, preK, j, &prp);
    for (; j < end; ++j)
      __compute1(K, preK, j, &prp);
  }

  template<typename FP>
  void compute(FP *K, const FP* preK) {
    std::vector<std::future<void>> fut;
    int width = n / threads;
    for (int i = 0; i < threads - 1; ++i) {
      int start = i * width;
      int end = std::min((i + 1) * width, n);
      fut.push_back(pool->enqueue([this, K, preK, start, end]() { task(K, preK, start, end); }));
    }
    int start = (threads-1) * width;
    int end = std::min((threads + 1) * width, n);
    task(K, preK, start, end);

    for (auto &f : fut)
      f.get();
  }

};
#endif
