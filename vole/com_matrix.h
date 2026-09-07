#ifndef DIST_PSI_COM_MATRIX_FP_H__
#define DIST_PSI_COM_MATRIX_FP_H__

// Commitment matrix for the committed VOLE (paper Sec. 4.1, dual-LPN commitment)
//
//   com = H_u * e_u + H_r * e_r,   H_u in F^{n_com x N},  H_r in F^{n_com x N_com}
//
// applied to the RAW regular-sparse vector [e_u || e_r] (no accumulator on the
// commitment path). Both blocks are random matrices regenerated on the fly from
// a PRP, column by column:
//
//   * H_u : column-sparse random (col_weight 32 by default). Column j has `col_weight` entries at rows
//           derived from PRP_u(j), each with a random coefficient. Every column
//           is hit (a two-position move of a nonzero always changes com), and a
//           support of size <= 2(t + t_com) is generically full column rank
//           (no second opening on the same or a different support). With
//           dense = true the column is fully dense (Theorem 1 of the paper as
//           written, ~10x the cost).
//   * H_r : dense random (the paper's H_2), so H_r * e_r is a dual-LPN sample
//           that masks H_u * e_u and gives computational hiding.
//
// Column generation: one PRP block per (row, coefficient) pair. The row index is
// taken from the high 32 bits of the block and the coefficient from
// FP::from_block on the whole block, so for 128-bit fields the coefficient
// shares bits with the row selector; conditioned on the row it still has >= 96
// uniform bits, which is all the rank argument needs.
//
// The same linear map is applied to the verifier's key vector, the committer's
// value||MAC vector (T = FP or FPS), and -- via compute_sparse -- to the king's
// sparse value vector.

#include "emp-tool/emp-tool.h"
#include "vole/fields/field_config.h"
#include <vector>

template <typename FP>
class ComMatrixFp {
public:
  std::size_t n_com, N, N_com, col_weight;
  bool dense;
  emp::block seed_u, seed_r;
  ThreadPool *pool;
  std::size_t threads;

  ComMatrixFp(std::size_t n_com, std::size_t N, std::size_t N_com,
              std::size_t col_weight, emp::block seed_u, emp::block seed_r,
              ThreadPool *pool, std::size_t threads, bool dense = false)
      : n_com(n_com), N(N), N_com(N_com), col_weight(col_weight), dense(dense),
        seed_u(seed_u), seed_r(seed_r), pool(pool), threads(threads) {
    if (col_weight > n_com) this->col_weight = n_com;
    if (threads < 1) this->threads = 1;
  }

  // Number of (row, coef) entries in column j.
  std::size_t width(std::size_t j) const {
    return (j >= N || dense) ? n_com : col_weight;
  }
  std::size_t max_width() const { return n_com; }

  // Per-thread scratch: a PRP for each block and a block buffer.
  struct Ctx {
    emp::PRP prp_u, prp_r;
    std::vector<emp::block> buf;
    Ctx(emp::block su, emp::block sr, std::size_t w) : prp_u(su), prp_r(sr), buf(w) {}
  };

  // Materialize column j: rows[k], coef[k] for k < width(j). Rows may repeat
  // (a repeat just merges two coefficients); coefficients are uniform in F.
  std::size_t column(std::size_t j, Ctx &ctx, uint32_t *rows, FP *coef) const {
    std::size_t w = width(j);
    emp::block *b = ctx.buf.data();
    for (std::size_t m = 0; m < w; ++m) b[m] = emp::makeBlock((uint64_t)j, (uint64_t)m);
    if (j < N) ctx.prp_u.permute_block(b, (int)w);
    else       ctx.prp_r.permute_block(b, (int)w);
    bool full = (w == n_com);
    for (std::size_t m = 0; m < w; ++m) {
      if (full) rows[m] = (uint32_t)m;
      else rows[m] = (uint32_t)(((uint64_t)_mm_extract_epi64(b[m], 1) >> 32) % n_com);
      coef[m].from_block(b[m]);
    }
    return w;
  }

  // out[0..n_com) = H * in   (in has length N + N_com; T = FP or FPS)
  template <typename T>
  void compute(T *out, const T *in) const {
    std::size_t total = N + N_com;
    std::size_t nth = std::min(threads, total);
    std::vector<std::vector<T>> acc(nth, std::vector<T>(n_com));
    std::vector<std::future<void>> fut;
    // Partition columns by cost (entries per column), not by count: the dense
    // H_r columns at the end cost n_com each versus col_weight for H_u, and
    // with the ring's large n_com they dominate.
    const double cost_u = (double)width(0), cost_r = (double)n_com;
    const double total_cost = (double)N * cost_u + (double)N_com * cost_r;
    auto bound = [&](std::size_t th) -> std::size_t {
      if (th >= nth) return total;
      double target = total_cost * (double)th / (double)nth;
      double cu = (double)N * cost_u;
      if (target <= cu) return std::min(N, (std::size_t)(target / cost_u));
      return std::min(total, N + (std::size_t)((target - cu) / cost_r));
    };
    for (std::size_t th = 0; th < nth; ++th) {
      std::size_t s = bound(th), e = bound(th + 1);
      auto job = [this, th, s, e, in, &acc]() {
        Ctx ctx(seed_u, seed_r, max_width());
        std::vector<uint32_t> rows(max_width());
        std::vector<FP> coef(max_width());
        T *a = acc[th].data();
        for (std::size_t i = 0; i < n_com; ++i) a[i].setZero();
        for (std::size_t j = s; j < e; ++j) {
          std::size_t w = column(j, ctx, rows.data(), coef.data());
          const T &v = in[j];
          for (std::size_t k = 0; k < w; ++k)
            a[rows[k]] = a[rows[k]] + v * coef[k];
        }
      };
      if (th + 1 < nth) fut.push_back(pool->enqueue(job));
      else job();
    }
    for (auto &f : fut) f.get();
    for (std::size_t i = 0; i < n_com; ++i) {
      T s = acc[0][i];
      for (std::size_t th = 1; th < nth; ++th) s = s + acc[th][i];
      out[i] = s;
    }
  }

  // out = H * (sparse vector with val[k] at position pos[k]); positions index the
  // concatenation [e_u || e_r]. Cost cnt * col_weight, used by the king.
  void compute_sparse(FP *out, const std::size_t *pos, const FP *val,
                      std::size_t cnt) const {
    Ctx ctx(seed_u, seed_r, max_width());
    std::vector<uint32_t> rows(max_width());
    std::vector<FP> coef(max_width());
    for (std::size_t i = 0; i < n_com; ++i) out[i].setZero();
    for (std::size_t k = 0; k < cnt; ++k) {
      std::size_t w = column(pos[k], ctx, rows.data(), coef.data());
      for (std::size_t m = 0; m < w; ++m)
        out[rows[m]] = out[rows[m]] + val[k] * coef[m];
    }
  }
};

#endif  // DIST_PSI_COM_MATRIX_FP_H__
