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

    k_mask = 1;
    while (k_mask < (uint32_t)k) {
      k_mask <<= 1;
      k_mask = k_mask | 0x1;
    }
  }

  template<typename FP>
  void __compute4(FP *K, const FP* preK, int i, PRP *prp) {
    block tmp[10];
    for (int m = 0; m < 10; ++m)
      tmp[m] = makeBlock(i, m);
    prp->permute_block(tmp, 10);
    int *index = (int *)(tmp);
    for (int j = 0; j < 40; ++j) {
      index[j] = index[j] & k_mask;
      index[j] = (index[j] >= k) ? index[j] - k : index[j];
    }

    // Lazy reduction: raw adds, reduced every FP::lazy_adds terms (5 for the
    // 61-bit Mersenne field, i.e. 2 reductions per output instead of 10).
    FP tmpv[4];
    tmpv[0] = K[i];
    tmpv[1] = K[i+1];
    tmpv[2] = K[i+2];
    tmpv[3] = K[i+3];
    int *p = (int *)(tmp);
    unsigned pending = 0;
    for (int j = 0; j < 10; ++j) {
      tmpv[0].add_raw(preK[*(p++)]);
      tmpv[1].add_raw(preK[*(p++)]);
      tmpv[2].add_raw(preK[*(p++)]);
      tmpv[3].add_raw(preK[*(p++)]);
      if (++pending == FP::lazy_adds) {
        tmpv[0].reduce(); tmpv[1].reduce(); tmpv[2].reduce(); tmpv[3].reduce();
        pending = 0;
      }
    }
    if (pending) { tmpv[0].reduce(); tmpv[1].reduce(); tmpv[2].reduce(); tmpv[3].reduce(); }
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
    prp->permute_block(tmp, 3);
    int *r = (int *)(tmp);
    for (int j = 0; j < 10; ++j) {
      r[j] = r[j] & k_mask;
      r[j] = (r[j] >= k) ? r[j] - k : r[j];
    }

    FP tmpv = K[i];
    unsigned pending = 0;
    for(int j = 0; j < 10; ++j) {
      tmpv.add_raw(preK[r[j]]);
      if (++pending == FP::lazy_adds) { tmpv.reduce(); pending = 0; }
    }
    if (pending) tmpv.reduce();
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
    vector<std::future<void>> fut;
    int width = n / threads;
    for (int i = 0; i < threads - 1; ++i) {
      int start = i * width;
      int end = min((i + 1) * width, n);
      fut.push_back(pool->enqueue([this, K, preK, start, end]() { task(K, preK, start, end); }));
    }
    int start = (threads-1) * width;
    int end = min((threads + 1) * width, n);
    task(K, preK, start, end);

    for (auto &f : fut)
      f.get();
  }

};
#endif
