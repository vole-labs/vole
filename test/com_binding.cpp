// Local (non-networked) checks on the committed-VOLE commitment layer,
// com = H_u * e_u + H_r * e_r (vole/com_matrix.h), using the king-local mode of
// CVoleFp to obtain an honest committer's sparse vector, x and com.
//
//  1) hiding sanity : no com_i equals x_i, and no com_i - x_i is a single +-beta.
//  2) shift-by-one  : moving one nonzero to a neighbouring position (same
//                     value) must always change com. (With the old EA-code
//                     commitment this succeeded 99% of the time at 2^20.)
//  3) rank          : the columns of [H_u | H_r] on the honest support, and on
//                     the union of the honest support with a second random
//                     regular support (2(t + t_com) columns), are linearly
//                     independent, so no second opening exists on any such
//                     support (paper Thm. 1). Fields only (needs inverses).
//
// usage: test_com_binding            (runs fp61, fp107, f2k, z2k)
#include "emp-tool/emp-tool.h"
#include "vole/cvole.h"
#include "vole/fields/fp61.h"
#include "vole/fields/fp61x2.h"
#include "vole/fields/fp107.h"
#include "vole/fields/fp107x2.h"
#include "vole/fields/fp2x128.h"
#include "vole/fields/z2k.h"
#include <vector>
#include <string>
#include <type_traits>

using namespace emp;

static int failures = 0;
#define CHECK(cond, msg)                                                   \
  do {                                                                     \
    if (!(cond)) { printf("  FAIL: %s\n", msg); ++failures; }              \
    else         { printf("  ok  : %s\n", msg); }                          \
  } while (0)

// Column rank of the n_com x |cols| matrix whose columns are given, over a
// field (uses FP::inv). Reduces each new column against the pivots in order.
template <typename FP>
std::size_t column_rank(std::vector<std::vector<FP>> cols, std::size_t n) {
  std::vector<std::size_t> prow;
  std::vector<std::vector<FP>> piv;
  for (auto &v : cols) {
    for (std::size_t k = 0; k < piv.size(); ++k) {
      FP c = v[prow[k]];
      if (c == FP(0)) continue;
      for (std::size_t i = 0; i < n; ++i) v[i] = v[i] - piv[k][i] * c;
    }
    std::size_t r = n;
    for (std::size_t i = 0; i < n; ++i) if (!(v[i] == FP(0))) { r = i; break; }
    if (r == n) continue;  // dependent
    FP inv = v[r].inv();
    for (std::size_t i = 0; i < n; ++i) v[i] = v[i] * inv;
    prow.push_back(r);
    piv.push_back(v);
  }
  return piv.size();
}

template <typename FP>
std::vector<FP> column_of(const ComMatrixFp<FP> &H, std::size_t j) {
  typename ComMatrixFp<FP>::Ctx ctx(H.seed_u, H.seed_r, H.max_width());
  std::vector<uint32_t> rows(H.max_width());
  std::vector<FP> coef(H.max_width());
  std::size_t w = H.column(j, ctx, rows.data(), coef.data());
  std::vector<FP> v(H.n_com);
  for (std::size_t k = 0; k < w; ++k) v[rows[k]] = v[rows[k]] + coef[k];
  return v;
}

// Rank checks (fields only). Dispatched on field_traits so the ring type never
// instantiates FP::inv().
template <typename FP, typename FPS>
void rank_checks(CVoleFp<NetIO, FP, FPS> &king, const CVoleFpParam &P,
                 const std::vector<std::size_t> &pos, std::true_type) {
  std::size_t N = P.N(), leave_u = 1ull << P.log_bin, leave_r = 1ull << P.log_bin_com;
  std::size_t s = pos.size();
  std::vector<std::vector<FP>> cols;
  for (std::size_t j = 0; j < s; ++j) cols.push_back(column_of(*king.hcom, pos[j]));
  std::size_t r1 = column_rank(cols, P.n_com);
  printf("  rank on honest support: %zu / %zu\n", r1, s);
  CHECK(r1 == s, "commitment map injective on the honest support");

  PRG prg;
  for (int trial = 0; trial < 3; ++trial) {
    std::vector<std::vector<FP>> cols2 = cols;
    std::vector<std::size_t> p2;
    for (std::size_t j = 0; j < s; ++j) {
      std::size_t base = (j < P.t) ? j * leave_u : N + (j - P.t) * leave_r;
      std::size_t leave = (j < P.t) ? leave_u : leave_r;
      uint64_t r; prg.random_data(&r, 8);
      std::size_t q = base + (r % leave);
      if (q == pos[j]) q = base + ((r + 1) % leave);
      p2.push_back(q);
    }
    for (auto q : p2) cols2.push_back(column_of(*king.hcom, q));
    std::size_t r2 = column_rank(cols2, P.n_com);
    printf("  rank on double support (trial %d): %zu / %zu\n", trial, r2, 2 * s);
    CHECK(r2 == 2 * s, "no second opening on a double support");
  }
}
template <typename FP, typename FPS>
void rank_checks(CVoleFp<NetIO, FP, FPS> &, const CVoleFpParam &,
                 const std::vector<std::size_t> &, std::false_type) {
  printf("  (ring type: rank test skipped, no inverses)\n");
}

template <typename FP, typename FPS>
void run(const char *tag) {
  printf("[%s]\n", tag);
  CVoleFpParam P = cvole_resolve_param<FP>(cvole_fp_n2to14);
  std::size_t N = P.N(), leave_u = 1ull << P.log_bin, leave_r = 1ull << P.log_bin_com;
  printf("  n_com (derived for this field) = %zu\n", P.n_com);

  CVoleFp<NetIO, FP, FPS> king(4, P);
  king.setup_prog(makeBlock(0xC0FFEEULL, 0x1234ULL));
  king.setup_local();
  std::vector<FP> x(P.n), com(P.n_com);
  king.extend_local(x.data(), com.data());
  std::vector<std::size_t> pos = king.dbg_pos;
  std::vector<FP> val = king.dbg_val;
  std::size_t s = pos.size();  // t + t_com
  CHECK(s == P.t + P.t_com, "honest support has t + t_com nonzeros");

  // 1) hiding sanity
  std::size_t eq = 0, single = 0;
  for (std::size_t i = 0; i < std::min(P.n, P.n_com); ++i) {
    if (com[i] == x[i]) { ++eq; continue; }
    FP d = com[i] - x[i];
    for (std::size_t j = 0; j < s; ++j)
      if (d == val[j] || d == val[j].negate()) { ++single; break; }
  }
  CHECK(eq == 0, "no com_i equals x_i");
  CHECK(single == 0, "no com_i - x_i equals a single +-beta");

  // 2) shift-by-one, every block of e_u and e_r
  std::size_t changed = 0;
  std::vector<FP> com2(P.n_com);
  for (std::size_t j = 0; j < s; ++j) {
    std::vector<std::size_t> p2 = pos;
    std::size_t base = (j < P.t) ? j * leave_u : N + (j - P.t) * leave_r;
    std::size_t leave = (j < P.t) ? leave_u : leave_r;
    std::size_t off = pos[j] - base;
    p2[j] = base + (off + 1 < leave ? off + 1 : off - 1);
    king.hcom->compute_sparse(com2.data(), p2.data(), val.data(), s);
    bool same = true;
    for (std::size_t i = 0; i < P.n_com; ++i) if (!(com2[i] == com[i])) same = false;
    if (!same) ++changed;
  }
  printf("  shift-by-one changed com in %zu / %zu blocks\n", changed, s);
  CHECK(changed == s, "shift-by-one always changes com");

  // 3) rank of the commitment map on the honest / double support (fields)
  rank_checks(king, P, pos, std::integral_constant<bool, !field_traits<FP>::is_ring>{});
}

int main(int argc, char **argv) {
  run<FP61, FP61x2>("fp61");
  run<FP107, FP107x2>("fp107");
  run<FP2x128, FP2x128x2>("f2k");
  run<Z2k64, Z2k64x2>("z2k");
  if (failures) { printf("%d check(s) FAILED\n", failures); return 1; }
  printf("all commitment checks passed\n");
  return 0;
}
