// Local (no network) differential checks of the MPFSS check arithmetic
// (vole/fields/utils.h) against the coefficient-array form it replaced:
//  * field_poly_hash(seed, n, x) == <field_uni_hash_coeff_gen(seed, n), x>
//    for every field, lengths around the 8-way chain split and the real tree
//    width, and the edge seeds 0, 1, -1 (where the field has them);
//  * field_hash_coeff_fixed(seed, e, bits) == chi[e] for every e < 2^bits
//    (bits = 6) and random, 0 and maximal e at bits = 12, 13.
//
// usage: test_poly_hash
#include "emp-tool/emp-tool.h"
#include "vole/fields/fp61.h"
#include "vole/fields/fp107.h"
#include "vole/fields/fp192.h"
#include "vole/fields/fp25519.h"
#include "vole/fields/fpp256.h"
#include "vole/fields/fp2x128.h"
#include "vole/fields/f2.h"
#include <vector>

using namespace emp;

static int failures = 0;
#define CHECK(cond, msg)                                                   \
  do {                                                                     \
    if (!(cond)) { printf("  FAIL: %s\n", msg); ++failures; }              \
    else         { printf("  ok  : %s\n", msg); }                          \
  } while (0)

template <typename FP>
static FP old_hash(FP seed, const std::vector<FP> &x) {
  std::vector<FP> chi(x.size());
  field_uni_hash_coeff_gen(chi.data(), seed, (int)x.size());
  return field_inn_prdt_sum_red(chi.data(), x.data(), (int)x.size());
}

template <typename FP>
static void run(const char *tag, std::vector<FP> edge_seeds) {
  printf("[%s]\n", tag);
  PRG prg;
  auto rnd = [&]() {
    FP v;
    if constexpr (requires { v.rand(prg); }) v.rand(prg);
    else { block b; prg.random_block(&b, 1); v.from_block(b); }
    return v;
  };
  for (int i = 0; i < 6; ++i) edge_seeds.push_back(rnd());

  const std::size_t lens[] = {1, 2, 3, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65, 100, 512, 4096};
  bool ok_hash = true;
  for (const FP &seed : edge_seeds)
    for (std::size_t n : lens) {
      std::vector<FP> x(n);
      for (auto &e : x) e = rnd();
      if (n % 3 == 0) x.back().setZero();  // zero top term
      FP got = field_poly_hash(seed, n, [&](std::size_t i) { return x[i]; });
      if (got != old_hash(seed, x)) ok_hash = false;
    }
  CHECK(ok_hash, "Horner form == coefficient array + inner product (all lengths, edge seeds)");

  bool ok_pow = true;
  for (const FP &seed : edge_seeds) {
    const int bits = 6;
    std::vector<FP> chi(1 << bits);
    field_uni_hash_coeff_gen(chi.data(), seed, 1 << bits);
    for (std::size_t e = 0; e < chi.size(); ++e)
      if (field_hash_coeff_fixed(seed, e, bits) != chi[e]) ok_pow = false;
    for (int big : {12, 13}) {
      std::vector<FP> c(std::size_t{1} << big);
      field_uni_hash_coeff_gen(c.data(), seed, 1 << big);
      for (int t = 0; t < 64; ++t) {
        uint32_t e;
        prg.random_data_unaligned(&e, sizeof(e));
        e &= (1u << big) - 1;
        if (t == 0) e = 0;
        if (t == 1) e = (1u << big) - 1;
        if (field_hash_coeff_fixed(seed, e, big) != c[e]) ok_pow = false;
      }
    }
  }
  CHECK(ok_pow, "fixed-sequence chi[e] == array chi[e] (all e < 2^6; random, 0, max e at 2^12, 2^13)");
}

int main() {
  run<FP61>("FP61", {FP61((uint64_t)0), FP61((uint64_t)1), FP61((uint64_t)FP61::PR - 1)});
  run<FP107>("FP107", {FP107((uint64_t)0), FP107((uint64_t)1), FP107((uint64_t)0) - FP107((uint64_t)1)});
  run<FP192>("FP192", {FP192((uint64_t)0), FP192((uint64_t)1), FP192(-1)});
  run<FP25519>("FP25519", {FP25519((uint64_t)0), FP25519((uint64_t)1), FP25519(-1)});
  run<FPP256>("FPP256", {FPP256((uint64_t)0), FPP256((uint64_t)1), FPP256(-1)});
  run<FP2x128>("FP2x128", {FP2x128((uint64_t)0), FP2x128((uint64_t)1)});
  run<F2kKey>("F2kKey", {F2kKey((uint64_t)0), F2kKey((uint64_t)2)});
  if (failures) { printf("%d check(s) FAILED\n", failures); return 1; }
  printf("all poly-hash checks passed\n");
  return 0;
}
