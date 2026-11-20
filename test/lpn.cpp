// Local (no network) checks of the primal-LPN expansion (vole/lpn.h) and the
// parameter sets (vole/vole_triple.h):
//  * the five sets are emp-ot's ferret_b9 ... ferret_b13 (t, k = 2^logk,
//    2^tree_depth, N), with two ferret_b9 bootstrap stages that cover each
//    set's t + k + 1 base pairs;
//  * the public matrix depends on the root, the stage domain and the round,
//    and on nothing else (not the thread count);
//  * the expansion matches a one-entry-at-a-time reference of the stream
//    layout: unit entries for prime fields and the F_2-valued key type,
//    one pseudorandom GF(2^128) coefficient per entry for FP2x128;
//  * over GF(2^128) one bit plane of the secret reaches other bit planes of
//    the output (no 128-way split), while F2kKey stays on its plane;
//  * MAC = key + value * Delta survives the expansion (FP61, GF(2^128)).
//
// usage: test_lpn
#include "emp-tool/emp-tool.h"
#include "vole/vole-fp.h"
#include "vole/fields/field_config.h"
#include <vector>

using namespace emp;

static int failures = 0;
#define CHECK(cond, msg)                                                   \
  do {                                                                     \
    if (!(cond)) { printf("  FAIL: %s\n", msg); ++failures; }              \
    else         { printf("  ok  : %s\n", msg); }                          \
  } while (0)

// ---- parameter sets ----------------------------------------------------
static void check_sets() {
  printf("[parameter sets]\n");
  struct Row { const char *name; const PrimalLPNParameterFp61 *p; std::size_t t, logk, depth, N; };
  const Row rows[] = {
      {"ferret_b9", &ferret_b9, 2400, 16, 9, 1228800},
      {"ferret_b10", &ferret_b10, 2800, 17, 10, 2867200},
      {"ferret_b11", &ferret_b11, 3800, 18, 11, 7782400},
      {"ferret_b12", &ferret_b12, 4200, 19, 12, 17203200},
      {"ferret_b13", &ferret_b13, 4500, 20, 13, 36864000},
  };
  bool ok = true;
  for (const Row &r : rows) {
    const PrimalLPNParameterFp61 &p = *r.p;
    bool m = p.t == r.t && p.k == (std::size_t{1} << r.logk) && p.log_bin_sz == r.depth && p.n == r.N &&
             p.t_pre == 2400 && p.k_pre == 65536 && p.log_bin_sz_pre == 9 && p.n_pre == 1228800 &&
             p.t_pre0 == 2400 && p.k_pre0 == 65536 && p.log_bin_sz_pre0 == 9 && p.n_pre0 == 1228800 &&
             p.n_pre >= p.k + p.t + 1 && p.n_pre0 >= p.k_pre + p.t_pre + 1;
    if (!m) printf("    mismatch in %s\n", r.name);
    ok = ok && m;
  }
  CHECK(ok, "b9..b13 = (t, logk, tree_depth) of emp-ot's sets; bootstrap stages cover t + k + 1");
  CHECK(fp_default.n == ferret_b12.n && fp_default.k == ferret_b12.k && fp_default.t == ferret_b12.t,
        "default set is ferret_b12");
  CHECK(PrimalLPNParameterFp61::is_pow2(65536) && !PrimalLPNParameterFp61::is_pow2(158000),
        "power-of-two check accepts 2^16, rejects 158000");
}

// ---- reference expansion ------------------------------------------------
// One entry at a time, straight from the stream definition in lpn.h.
template <typename FP>
static void reference(std::vector<FP> &out, const std::vector<FP> &pre, int n, int k, block root,
                      const char *domain, uint64_t round) {
  const block ik = RO(domain, root).absorb(std::string_view("index")).absorb(round).squeeze_block();
  const block ck = RO(domain, root).absorb(std::string_view("coefficient")).absorb(round).squeeze_block();
  PRP prp(ik), cprp(ck);
  for (int o = 0; o < n; ++o) {
    const int b = o / 32, pos = o % 32;
    block idxb[80];
    for (int m = 0; m < 80; ++m) idxb[m] = makeBlock((uint64_t)b, (uint64_t)m);
    prp.permute_block(idxb, 80);
    const uint32_t *w = (const uint32_t *)idxb;
    for (int j = 0; j < 10; ++j) {
      const FP &e = pre[w[pos * 10 + j] & (uint32_t)(k - 1)];
      if constexpr (lpn_coefficients<FP>::value) {
        block c = makeBlock((uint64_t)b, (uint64_t)(pos * 10 + j));
        cprp.permute_block(&c, 1);
        lpn_coefficient_t<FP> a;
        a.from_block(c);
        out[o] = out[o] + e * a;
      } else {
        out[o] = out[o] + e;
      }
    }
  }
}

template <typename FP>
static std::vector<FP> rand_vec(PRG &prg, int len) {
  std::vector<FP> v(len);
  for (auto &e : v) {
    if constexpr (requires { e.rand(prg); }) e.rand(prg);
    else { FP2x128 lo, hi; lo.rand(prg); hi.rand(prg); e.setLowHigh(lo, hi); }  // FP2x128x2
  }
  return v;
}

template <typename FP>
static std::vector<FP> expand(ThreadPool *pool, int threads, int n, int k, block root, const char *domain,
                              const std::vector<FP> &init, const std::vector<FP> &pre, int rounds = 1) {
  LpnFp<10> lpn(n, k, pool, threads, root, domain);
  std::vector<FP> out = init;
  for (int r = 0; r < rounds; ++r) { out = init; lpn.compute(out.data(), pre.data()); }
  return out;  // the output of round `rounds - 1`
}

template <typename FP>
static void check_matrix(ThreadPool *pool, const char *tag) {
  printf("[%s]\n", tag);
  PRG prg;
  const int n = 4096 + 17, k = 1024;  // a partial last batch
  block root, root2;
  prg.random_block(&root, 1);
  prg.random_block(&root2, 1);
  std::vector<FP> init = rand_vec<FP>(prg, n), pre = rand_vec<FP>(prg, k);

  std::vector<FP> a = expand(pool, 4, n, k, root, kLpnMain, init, pre);
  std::vector<FP> ref = init;
  reference(ref, pre, n, k, root, kLpnMain, 0);
  CHECK(a == ref, "matches the one-entry-at-a-time reference (round 0)");

  std::vector<FP> a1 = expand(pool, 1, n, k, root, kLpnMain, init, pre);
  CHECK(a == a1, "independent of the thread count (4 vs 1)");

  std::vector<FP> r1 = expand(pool, 4, n, k, root, kLpnMain, init, pre, 2);
  std::vector<FP> ref1 = init;
  reference(ref1, pre, n, k, root, kLpnMain, 1);
  CHECK(r1 == ref1 && r1 != a, "round 1 reads a fresh matrix (and matches the reference)");

  CHECK(expand(pool, 4, n, k, root2, kLpnMain, init, pre) != a, "a different root gives a different matrix");
  CHECK(expand(pool, 4, n, k, root, kLpnPre, init, pre) != a &&
        expand(pool, 4, n, k, root, kLpnPre0, init, pre) != a, "each stage domain gives a different matrix");
}

// One bit plane of the secret: pre[j] = 1 for every j, so a unit-entry output
// is a 0/1 sum and stays on bit plane 0; a coefficient output does not.
template <typename FP>
static bool leaves_plane0(ThreadPool *pool) {
  const int n = 1024, k = 256;
  std::vector<FP> init(n), pre(k);
  for (auto &e : init) e.setZero();
  for (auto &e : pre) e = FP((uint64_t)1);
  block root = makeBlock(1, 2);
  std::vector<FP> out = expand(pool, 2, n, k, root, kLpnMain, init, pre);
  for (auto &e : out)
    if (e != FP((uint64_t)0) && e != FP((uint64_t)1)) return true;
  return false;
}

// MAC = key + value * Delta after expanding keys (FP) and bundles (FPS).
template <typename FP, typename FPS>
static bool mac_relation(ThreadPool *pool) {
  PRG prg;
  const int n = 2048 + 5, k = 512;
  FP delta;
  delta.rand(prg);
  block root;
  prg.random_block(&root, 1);
  std::vector<FP> key0 = rand_vec<FP>(prg, n), keyp = rand_vec<FP>(prg, k);
  std::vector<FP> val0 = rand_vec<FP>(prg, n), valp = rand_vec<FP>(prg, k);
  std::vector<FPS> b0(n), bp(k);
  for (int i = 0; i < n; ++i) b0[i].setLowHigh(key0[i] + val0[i] * delta, val0[i]);
  for (int i = 0; i < k; ++i) bp[i].setLowHigh(keyp[i] + valp[i] * delta, valp[i]);
  std::vector<FP> key = expand(pool, 3, n, k, root, kLpnMain, key0, keyp, 2);
  std::vector<FPS> bun = expand(pool, 3, n, k, root, kLpnMain, b0, bp, 2);
  for (int i = 0; i < n; ++i)
    if (bun[i].getLow() != key[i] + bun[i].getHigh() * delta) return false;
  return true;
}

int main() {
  ThreadPool pool(4);
  check_sets();
  check_matrix<FP61>(&pool, "FP61 (unit entries)");
  check_matrix<FP192>(&pool, "FP192 (unit entries, lazy accumulator)");
  check_matrix<FP25519>(&pool, "FP25519 (unit entries, lazy accumulator)");
  check_matrix<FPP256>(&pool, "FPP256 (unit entries, lazy accumulator)");
  check_matrix<F2kKey>(&pool, "F2kKey (F_2 values, unit entries)");
  check_matrix<FP2x128>(&pool, "FP2x128 (GF(2^128) coefficients)");
  check_matrix<FP2x128x2>(&pool, "FP2x128x2 (GF(2^128) coefficients, value/MAC lanes)");

  printf("[bit planes / MAC relation]\n");
  CHECK(leaves_plane0<FP2x128>(&pool), "GF(2^128): one bit plane of the secret reaches other planes");
  CHECK(!leaves_plane0<F2kKey>(&pool), "F2kKey: unit entries keep F_2 values on plane 0");
  CHECK((mac_relation<FP2x128, FP2x128x2>(&pool)), "GF(2^128): MAC = key + value*Delta after two rounds");
  CHECK((mac_relation<FP61, FP61x2>(&pool)), "FP61: MAC = key + value*Delta after two rounds");
  CHECK((mac_relation<FPP256, FPP256x2>(&pool)), "P-256: MAC = key + value*Delta after two rounds");

  if (failures) { printf("%d check(s) FAILED\n", failures); return 1; }
  printf("all LPN checks passed\n");
  return 0;
}
