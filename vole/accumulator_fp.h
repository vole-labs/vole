#ifndef DIST_PSI_ACCUMULATOR_FP_H__
#define DIST_PSI_ACCUMULATOR_FP_H__

// Field-generic accumulator: the "accumulate" half of the expand-accumulate
// (EA) code [BCG+22], with field addition in place of XOR and templated on the
// element type (FP or FPS). Computes an in-place prefix sum:
// data[i] <- sum_{j<=i} data[j]. The two-pass version is the parallel layout;
// small/degenerate sizes fall back to a
// sequential pass (mathematically identical prefix sum).

#include "emp-tool/emp-tool.h"

class AccumFp {
public:
  std::size_t n;
  ThreadPool *pool;
  std::size_t threads;
  std::size_t width;

  AccumFp(std::size_t n, ThreadPool *pool, std::size_t threads)
      : n(n), pool(pool), threads(threads) {
    width = n / (threads + 2);
  }

  template <typename T>
  void first_pass(T *data) {
    std::vector<std::future<void>> fut;
    for (std::size_t i = 0; i < threads; ++i) {
      T *pointer = data + i * width;
      fut.push_back(pool->enqueue([this, pointer]() {
        for (std::size_t j = 1; j < width; ++j)
          pointer[j] = pointer[j - 1] + pointer[j];
      }));
    }
    T *pointer = data + threads * width;
    for (std::size_t j = 1; j < width; ++j)
      pointer[j] = pointer[j - 1] + pointer[j];

    for (auto &f : fut)
      f.get();
  }

  template <typename T>
  void second_pass(T *data) {
    std::vector<std::future<void>> fut;

    std::vector<T> block_sum(threads + 1);
    for (std::size_t i = 0; i < threads + 1; ++i) {
      block_sum[i] = data[(i + 1) * width - 1];
      if (i != 0)
        block_sum[i] = block_sum[i] + block_sum[i - 1];
    }

    for (std::size_t i = 0; i < threads; ++i) {
      T *pointer = data + (i + 1) * width;
      fut.push_back(pool->enqueue([this, i, pointer, &block_sum]() {
        for (std::size_t j = 0; j < width; ++j)
          pointer[j] = block_sum[i] + pointer[j];
      }));
    }
    std::size_t start = (threads + 1) * width;
    if (n > start) {
      T *pointer = data + start;
      pointer[0] = pointer[0] + block_sum[threads];
      for (std::size_t j = 1; j < n - start; ++j)
        pointer[j] = pointer[j - 1] + pointer[j];
    }

    for (auto &f : fut)
      f.get();
  }

  template <typename T>
  void compute(T *data) {
    if (width < 2 || threads == 0) {
      compute_single_pass(data);
      return;
    }
    first_pass(data);
    second_pass(data);
  }

  template <typename T>
  void compute_single_pass(T *data) {
    for (std::size_t i = 1; i < n; ++i)
      data[i] = data[i - 1] + data[i];
  }
};

#endif // DIST_PSI_ACCUMULATOR_FP_H__
