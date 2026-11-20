// Local (no network) checks of the 256-bit prime fields
// (vole/fields/fp25519.h: p = 2^255 - 19; vole/fields/fpp256.h: NIST P-256):
//  1) known answers: +, -, *, inv against Python big-int results
//     (fp256_kat.h), edge values and random pairs;
//  2) reduction: the U256 constructor on raw 256-bit values; P-256's 384-bit
//     sampling reduction (from_3blocks) on raw 384-bit values; 2^255-19's
//     from_2blocks (drop bit 255, reduce);
//  3) algebraic identities on random elements, rand() canonical;
//  4) the lazy accumulators Acc (FP and pair) against the reduced running sum;
//  5) pair lanes agree with the field.
//
// usage: test_fp256_arith
#include "emp-tool/emp-tool.h"
#include "vole/fields/fp25519.h"
#include "vole/fields/fp25519x2.h"
#include "vole/fields/fpp256.h"
#include "vole/fields/fpp256x2.h"
#include "test/fp256_kat.h"
#include <string>
#include <vector>

using namespace emp;

static int failures = 0;
#define CHECK(cond, msg)                                                   \
  do {                                                                     \
    if (!(cond)) { printf("  FAIL: %s\n", msg); ++failures; }              \
    else         { printf("  ok  : %s\n", msg); }                          \
  } while (0)

// "0x" + 16*nlimb hex digits -> little-endian limbs
static void parse(const char *hex, uint64_t *w, int nlimb) {
  std::string s(hex + 2);
  for (int l = 0; l < nlimb; ++l)
    w[l] = std::stoull(s.substr(s.size() - 16 * (l + 1), 16), nullptr, 16);
}
static U256 u256_of(const char *hex) { U256 r; parse(hex, r.w, 4); return r; }
template <typename F> static F raw(const U256 &v) { F r; r.val = v; return r; }  // no reduction

static bool canonical(const U256 &v, const U256 &p) {
  for (int i = 3; i >= 0; --i) {
    if (v.w[i] < p.w[i]) return true;
    if (v.w[i] > p.w[i]) return false;
  }
  return false;  // == p
}

template <typename F, std::size_t N>
static void known_answers(const char *(&kat)[N][6]) {
  std::size_t bad_add = 0, bad_sub = 0, bad_mul = 0, bad_inv = 0, bad_neg = 0;
  for (std::size_t i = 0; i < N; ++i) {
    F a = raw<F>(u256_of(kat[i][0])), b = raw<F>(u256_of(kat[i][1]));
    if (a + b != raw<F>(u256_of(kat[i][2]))) ++bad_add;
    if (a - b != raw<F>(u256_of(kat[i][3]))) ++bad_sub;
    if (a * b != raw<F>(u256_of(kat[i][4]))) ++bad_mul;
    if (a.inv() != raw<F>(u256_of(kat[i][5]))) ++bad_inv;
    if (b - a != -(a - b)) ++bad_neg;
  }
  printf("  %zu vectors\n", N);
  CHECK(bad_add == 0, "a + b");
  CHECK(bad_sub == 0, "a - b");
  CHECK(bad_mul == 0, "a * b");
  CHECK(bad_inv == 0, "inv(a) (0 -> 0)");
  CHECK(bad_neg == 0, "b - a == -(a - b)");
}

template <typename F, std::size_t N>
static void reductions256(const char *(&red)[N][2]) {
  std::size_t bad = 0;
  for (std::size_t i = 0; i < N; ++i)
    if (F(u256_of(red[i][0])) != raw<F>(u256_of(red[i][1]))) ++bad;
  printf("  %zu raw values\n", N);
  CHECK(bad == 0, "U256 constructor reduces any 256-bit value");
}

template <typename F, typename FS>
static void generic(const char *tag) {
  printf("[%s identities / sampling / accumulators]\n", tag);
  PRG prg;
  bool assoc = true, distr = true, inv = true, neg = true, shift = true, u64 = true, canon = true;
  for (int i = 0; i < 2000; ++i) {
    F a, b, c;
    a.rand(prg); b.rand(prg); c.rand(prg);
    if (!canonical(a.val, F::PR) || !canonical(b.val, F::PR)) canon = false;
    if ((a * b) * c != a * (b * c)) assoc = false;
    if (a * (b + c) != a * b + a * c) distr = false;
    if (!(a == F((uint64_t)0)) && a * a.inv() != F((uint64_t)1)) inv = false;
    if (a + (-a) != F((uint64_t)0) || a.negate() != -a) neg = false;
    int n = i % 256;
    F two_n((uint64_t)1);
    for (int k = 0; k < n; ++k) two_n = two_n + two_n;
    if ((a << n) != a * two_n) shift = false;
    uint64_t s = 0x9e3779b97f4a7c15ULL * (i + 1);
    if (a + s != a + F(s) || a * s != a * F(s) || a - s != a - F(s)) u64 = false;
  }
  CHECK(canon, "rand() always returns a canonical value");
  CHECK(assoc, "(ab)c = a(bc)");
  CHECK(distr, "a(b + c) = ab + ac");
  CHECK(inv, "a * inv(a) = 1");
  CHECK(neg, "a + (-a) = 0, negate() = -a");
  CHECK(shift, "a << n = a * 2^n for n < 256");
  CHECK(u64, "uint64 operand overloads");
  CHECK(F(-1) == F((uint64_t)0) - F((uint64_t)1) && F(-5) + F(5) == F((uint64_t)0), "int constructor (negatives)");

  // prg_elements: canonical, and the same stream gives the same elements
  {
    std::vector<F> x(4096), y(4096);
    block sd = makeBlock(7, 9);
    PRG g1(&sd), g2(&sd);
    F::prg_elements(g1, x.data(), x.size());
    F::prg_elements(g2, y.data(), y.size());
    bool ok = true;
    for (std::size_t i = 0; i < x.size(); ++i)
      if (!canonical(x[i].val, F::PR) || x[i] != y[i]) ok = false;
    CHECK(ok, "prg_elements: canonical and deterministic in the seed");
  }

  // Acc == reduced running sum, random and worst-case (p - 1) terms
  bool acc = true, acc2 = true;
  for (int trial = 0; trial < 40; ++trial) {
    int n = 1 + (trial * 127) % 5000;
    F start, run;
    if (trial % 2) start = F((uint64_t)0) - F((uint64_t)1); else start.rand(prg);
    run = start;
    typename F::Acc a(start);
    FS s0(start, start * F((uint64_t)3)), srun = s0;
    typename FS::Acc a2(s0);
    for (int i = 0; i < n; ++i) {
      F t;
      if (trial % 2) t = F((uint64_t)0) - F((uint64_t)1); else t.rand(prg);
      a.add(t); run = run + t;
      FS ts(t, t + t); a2.add(ts); srun = srun + ts;
    }
    if (a.get() != run) acc = false;
    if (a2.get() != srun) acc2 = false;
  }
  CHECK(acc, "F::Acc equals the reduced running sum (1..5000 terms, p-1 worst case)");
  CHECK(acc2, "pair Acc equals the reduced running sum");

  // pair lanes
  bool lanes = true;
  for (int i = 0; i < 500; ++i) {
    F a, b, c, d;
    a.rand(prg); b.rand(prg); c.rand(prg); d.rand(prg);
    FS x(a, b), y(c, d);
    if ((x + y).getLow() != a + c || (x + y).getHigh() != b + d) lanes = false;
    if ((x - y).getLow() != a - c || (x * y).getHigh() != b * d) lanes = false;
    if ((x * c).getLow() != a * c || x.negate().getHigh() != b.negate()) lanes = false;
    FS z = x; z.add_raw(y);
    if (z != x + y) lanes = false;
  }
  CHECK(lanes, "pair: lane-wise +, -, *, scalar ops, negate, add_raw");
}

int main() {
  printf("[2^255-19 known answers]\n");
  known_answers<FP25519>(KAT25519);
  printf("[2^255-19 reduction]\n");
  reductions256<FP25519>(RED25519);
  {
    // from_2blocks drops bit 255 and reduces: compare with the constructor
    PRG prg;
    bool ok = true, top = true;
    for (int i = 0; i < 20000; ++i) {
      block b[2];
      prg.random_block(b, 2);
      if (i == 0) { b[0] = makeBlock(~0ULL, ~0ULL); b[1] = makeBlock(~0ULL, ~0ULL); }
      if (i == 1) { b[0] = makeBlock(~0ULL, ~0ULL - 18); b[1] = makeBlock(~0ULL >> 1, ~0ULL); }  // = p
      FP25519 f; f.from_2blocks(b[0], b[1]);
      U256 v; std::memcpy(&v, b, 32); v.w[3] &= ~0ULL >> 1;
      if (f != FP25519(v)) ok = false;
      if (i > 1 && f.val.w[3] == 0) top = false;  // probability 2^-63 per draw
    }
    CHECK(ok, "from_2blocks = 255 low bits mod p (incl. all-ones and p itself)");
    CHECK(top, "from_2blocks uses the top limb");
  }
  generic<FP25519, FP25519x2>("2^255-19");

  printf("[P-256 known answers]\n");
  known_answers<FPP256>(KATP256);
  printf("[P-256 reduction]\n");
  reductions256<FPP256>(RED256);
  {
    std::size_t n = sizeof(RED384P256) / sizeof(RED384P256[0]), bad = 0;
    for (std::size_t i = 0; i < n; ++i) {
      uint64_t w[6];
      parse(RED384P256[i][0], w, 6);
      block b0, b1, b2;
      std::memcpy(&b0, &w[0], 16); std::memcpy(&b1, &w[2], 16); std::memcpy(&b2, &w[4], 16);
      FPP256 f; f.from_3blocks(b0, b1, b2);
      if (f != raw<FPP256>(u256_of(RED384P256[i][1]))) ++bad;
    }
    printf("  %zu raw 384-bit values\n", n);
    CHECK(bad == 0, "from_3blocks: 384-bit value mod p");
  }
  generic<FPP256, FPP256x2>("P-256");

  if (failures) { printf("%d check(s) FAILED\n", failures); return 1; }
  printf("all 256-bit field checks passed\n");
  return 0;
}
