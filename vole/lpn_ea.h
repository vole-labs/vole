#ifndef DIST_PSI_LPN_EA_FP_H__
#define DIST_PSI_LPN_EA_FP_H__

// Field-generic "expand" half of the expand-accumulate (EA) code [BCG+22]:
// each of the n outputs is the sum of `sparcity` entries of the length-k input
// chosen by a PRP, where sparcity = floor(ln(k) * 10) (expansion parameter
// C = 10, the conservative EA setting of the paper). Field addition replaces
// XOR, and `compute` is
// templated on the element type (FP or FPS) so the same matrix is applied to a
// key vector and to a packed value||MAC vector. Uses 64-bit indices (like
// LpnF2kSmall), so k may exceed 2^31.
//
// Semantics: compute(out, in) does out[i] += sum_j in[ idx_j(i) ]. Callers
// pre-zero `out` (or accumulate several expansions into the same `out`).

#include "emp-tool/emp-tool.h"
#include <cmath>

namespace emp {

class LpnFpEA {
public:
  int64_t k, n;
  int64_t sparcity;
  ThreadPool *pool;
  int threads;
  block seed;

  int64_t buf1_sz, buf4_sz;
  block **buf1 = nullptr;
  block **buf4 = nullptr;
  int64_t k_mask;

  LpnFpEA(int64_t n, int64_t k, ThreadPool *pool, int threads,
          block seed = zero_block) {
    this->n = n;
    this->k = k;
    this->pool = pool;
    this->threads = threads;
    this->seed = seed;

    // vole_f2k.h: sparcity = int64_t(log((double)k) * 10.0)  (natural log)
    sparcity = (int64_t)(std::log((double)k) * 10.0);
    if (sparcity < 1)
      sparcity = 1;

    buf1_sz = (sparcity + 1) / 2;
    buf4_sz = (sparcity * 4 + 1) / 2;
    buf1 = new block *[threads];
    for (int i = 0; i < threads; ++i)
      buf1[i] = new block[buf1_sz];
    buf4 = new block *[threads];
    for (int i = 0; i < threads; ++i)
      buf4[i] = new block[buf4_sz];

    k_mask = 1;
    while (k_mask < k) {
      k_mask <<= 1;
      k_mask = k_mask | 0x1;
    }
  }

  ~LpnFpEA() {
    for (int i = 0; i < threads; ++i)
      delete[] buf1[i];
    delete[] buf1;
    for (int i = 0; i < threads; ++i)
      delete[] buf4[i];
    delete[] buf4;
  }

  template <typename T>
  void add_one(T *out, const T *in, int64_t idx1, int64_t *idx2) {
    // Raw adds in groups of T::lazy_adds (compile-time constant), one reduce
    // per group; the remainder is reduced once at the end.
    T acc = out[idx1];
    const int L = (int)std::min<unsigned>(T::lazy_adds, 1u << 20);
    int j = 0;
    for (; j + L <= sparcity; j += L) {
      for (int jj = 0; jj < L; ++jj) acc.add_raw(in[idx2[j + jj]]);
      acc.reduce();
    }
    for (; j < sparcity; ++j) acc.add_raw(in[idx2[j]]);
    acc.reduce();
    out[idx1] = acc;
  }

  template <typename T>
  void __compute4(T *out, const T *in, int64_t i, PRP *prp, block *buf) {
    for (int m = 0; m < buf4_sz; ++m)
      buf[m] = makeBlock(i, m);
    prp->permute_block(buf, buf4_sz);
    int64_t *r = (int64_t *)(buf);
    for (int j = 0; j < 4 * sparcity; ++j) {
      r[j] = r[j] & k_mask;
      r[j] = r[j] >= k ? r[j] - k : r[j];
    }
    for (int64_t m = 0; m < 4; ++m)
      add_one(out, in, i + m, r + m * sparcity);
  }

  template <typename T>
  void __compute1(T *out, const T *in, int64_t i, PRP *prp, block *buf) {
    for (int m = 0; m < buf1_sz; ++m)
      buf[m] = makeBlock(i, m);
    prp->permute_block(buf, buf1_sz);
    int64_t *r = (int64_t *)(buf);
    for (int j = 0; j < sparcity; ++j) {
      r[j] = r[j] & k_mask;
      r[j] = r[j] >= k ? r[j] - k : r[j];
    }
    add_one(out, in, i, r);
  }

  template <typename T>
  void task(T *out, const T *in, int64_t start, int64_t end, block *b4,
            block *b1) {
    PRP prp(seed);
    int64_t j = start;
    for (; j < end - 4; j += 4)
      __compute4(out, in, j, &prp, b4);
    for (; j < end; ++j)
      __compute1(out, in, j, &prp, b1);
  }

  template <typename T>
  void compute(T *out, const T *in) {
    std::vector<std::future<void>> fut;
    int64_t width = (n + threads - 1) / threads;
    for (int i = 0; i < threads - 1; ++i) {
      int64_t start = i * width;
      int64_t end = std::min((i + 1) * width, n);
      if (start >= end)
        break;
      fut.push_back(pool->enqueue(
          [this, out, in, start, end, i]() {
            task(out, in, start, end, buf4[i], buf1[i]);
          }));
    }
    int64_t start = (int64_t)(threads - 1) * width;
    int64_t end = std::min((int64_t)threads * width, n);
    if (start < end)
      task(out, in, start, end, buf4[threads - 1], buf1[threads - 1]);

    for (auto &f : fut)
      f.get();
  }
};

} // namespace emp

#endif // DIST_PSI_LPN_EA_FP_H__
