// Local (no network) checks of the 128-bit prime field
// (vole/fields/fp128.h: p = 2^128 - 159):
//  1) known answers: +, -, *, inv against Python big-int results
//     (fp128_kat.h), edge values (0, 158, 159, p-1, 2^64, 2^127, ...) and
//     random pairs;
//  2) reduction: the U128 constructor and from_block on raw 128-bit values
//     (including [p, 2^128)); canonical() exactly on [0, p);
//  3) algebraic identities on random elements, rand() canonical, << n;
//  4) the lazy accumulators Acc (FP and pair) against the reduced running sum,
//     including long runs of p - 1;
//  5) pair lanes agree with the field.
//
// usage: test_fp128_arith
#include "emp-tool/emp-tool.h"
#include "vole/fields/fp128.h"
#include "vole/fields/fp128x2.h"
#include "test/fp128_kat.h"
#include <string>

using namespace emp;

static int failures = 0;
#define CHECK(cond, msg)                                                   \
  do {                                                                     \
    if (!(cond)) { printf("  FAIL: %s\n", msg); ++failures; }              \
    else         { printf("  ok  : %s\n", msg); }                          \
  } while (0)

static U128 u128_of(const char *hex) {
  std::string s(hex + 2);
  return U128{{std::stoull(s.substr(16, 16), nullptr, 16), std::stoull(s.substr(0, 16), nullptr, 16)}};
}
static FP128 raw(const U128 &v) { FP128 r; r.val = v; return r; }  // no reduction
static FP128 rnd(PRG &prg) { FP128 r; r.rand(prg); return r; }

static void known_answers() {
  printf("[known answers]\n");
  constexpr std::size_t N = sizeof(KAT128) / sizeof(KAT128[0]);
  std::size_t bad_add = 0, bad_sub = 0, bad_mul = 0, bad_inv = 0, bad_neg = 0, bad_can = 0;
  for (std::size_t i = 0; i < N; ++i) {
    FP128 a = raw(u128_of(KAT128[i][0])), b = raw(u128_of(KAT128[i][1]));
    FP128 s = a + b, d = a - b, m = a * b, v = a.inv();
    if (s != raw(u128_of(KAT128[i][2]))) ++bad_add;
    if (d != raw(u128_of(KAT128[i][3]))) ++bad_sub;
    if (m != raw(u128_of(KAT128[i][4]))) ++bad_mul;
    if (v != raw(u128_of(KAT128[i][5]))) ++bad_inv;
    if (b - a != -(a - b)) ++bad_neg;
    if (!s.canonical() || !d.canonical() || !m.canonical() || !v.canonical()) ++bad_can;
  }
  printf("  %zu vectors\n", N);
  CHECK(bad_add == 0, "a + b");
  CHECK(bad_sub == 0, "a - b");
  CHECK(bad_mul == 0, "a * b");
  CHECK(bad_inv == 0, "inv(a)");
  CHECK(bad_neg == 0, "b - a == -(a - b)");
  CHECK(bad_can == 0, "results canonical");
}

static void reductions() {
  printf("[reduction]\n");
  constexpr std::size_t N = sizeof(RED128) / sizeof(RED128[0]);
  std::size_t bad_ctor = 0, bad_blk = 0, bad_can = 0;
  for (std::size_t i = 0; i < N; ++i) {
    U128 in = u128_of(RED128[i][0]), want = u128_of(RED128[i][1]);
    if (FP128(in).val != want) ++bad_ctor;
    block b;
    std::memcpy(&b, &in, sizeof(b));
    FP128 f;
    f.from_block(b);
    if (f.val != want) ++bad_blk;
    if (raw(in).canonical() != (in == want)) ++bad_can;  // canonical iff in < p
  }
  printf("  %zu values\n", N);
  CHECK(bad_ctor == 0, "U128 constructor reduces mod p");
  CHECK(bad_blk == 0, "from_block reduces mod p");
  CHECK(bad_can == 0, "canonical() iff value < p");
  CHECK(raw(FP128::PR).canonical() == false && raw(U128{{~0ULL - 159, ~0ULL}}).canonical(),
        "canonical() boundary: p - 1 yes, p no");
}

static void identities() {
  printf("[identities]\n");
  PRG prg;
  bool ok_ring = true, ok_inv = true, ok_shift = true, ok_rand = true, ok_int = true;
  for (int i = 0; i < 20000; ++i) {
    FP128 a = rnd(prg), b = rnd(prg), c = rnd(prg);
    if (!a.canonical()) ok_rand = false;
    if ((a + b) * c != a * c + b * c) ok_ring = false;
    if (a * (b * c) != (a * b) * c) ok_ring = false;
    if (a - b + b != a || a + (-a) != FP128((uint64_t)0)) ok_ring = false;
    if (a != FP128((uint64_t)0) && a * a.inv() != FP128((uint64_t)1)) ok_inv = false;
    int n = i % 128;
    U128 two_n{{0, 0}};
    two_n.w[n >> 6] = 1ULL << (n & 63);
    if ((a << n) != a * raw(two_n)) ok_shift = false;
  }
  if (FP128(-1) != FP128((uint64_t)0) - FP128((uint64_t)1)) ok_int = false;
  if (FP128(-1).val != U128{{~0ULL - 159, ~0ULL}}) ok_int = false;
  CHECK(ok_rand, "rand() canonical");
  CHECK(ok_ring, "ring identities (distributive, associative, inverse of +)");
  CHECK(ok_inv, "a * inv(a) == 1");
  CHECK(ok_shift, "a << n == a * 2^n");
  CHECK(ok_int, "FP128(-1) == p - 1");
}

static void accumulators() {
  printf("[accumulators]\n");
  PRG prg;
  bool ok_rand = true, ok_max = true, ok_pair = true;
  for (int t = 0; t < 200; ++t) {
    FP128 x = rnd(prg), sum = x;
    FP128::Acc acc(x);
    for (int i = 0; i < 1 + t * 7; ++i) {
      FP128 y = rnd(prg);
      acc.add(y);
      sum = sum + y;
    }
    if (acc.get() != sum || !acc.get().canonical()) ok_rand = false;
  }
  const FP128 pm1 = FP128(-1);
  for (uint64_t cnt : {1ULL, 2ULL, 3ULL, 1000ULL, 1000000ULL}) {
    FP128::Acc acc(pm1);
    for (uint64_t i = 0; i < cnt; ++i) acc.add(pm1);
    // (cnt + 1) * (p - 1) = -(cnt + 1)
    if (acc.get() != -FP128(cnt + 1)) ok_max = false;
  }
  for (int t = 0; t < 100; ++t) {
    FP128x2 x(rnd(prg), rnd(prg)), sum = x;
    FP128x2::Acc acc(x);
    for (int i = 0; i < 50; ++i) {
      FP128x2 y(rnd(prg), rnd(prg));
      acc.add(y);
      sum = sum + y;
    }
    if (acc.get() != sum) ok_pair = false;
  }
  CHECK(ok_rand, "Acc == running sum (random terms)");
  CHECK(ok_max, "Acc on up to 10^6 + 1 copies of p - 1");
  CHECK(ok_pair, "FP128x2::Acc == running sum");
}

static void pairs() {
  printf("[pair lanes]\n");
  PRG prg;
  bool ok = true;
  for (int i = 0; i < 2000; ++i) {
    FP128 a0 = rnd(prg), a1 = rnd(prg), b0 = rnd(prg), b1 = rnd(prg), k = rnd(prg);
    FP128x2 a(a0, a1), b(b0, b1);
    if ((a + b) != FP128x2(a0 + b0, a1 + b1)) ok = false;
    if ((a - b) != FP128x2(a0 - b0, a1 - b1)) ok = false;
    if ((a * k) != FP128x2(a0 * k, a1 * k)) ok = false;
    if (a.negate() != FP128x2(-a0, -a1)) ok = false;
    if (a.getLow() != a0 || a.getHigh() != a1) ok = false;
  }
  CHECK(ok, "FP128x2 lane-wise +, -, * scalar, negate, low/high");
}

int main() {
  known_answers();
  reductions();
  identities();
  accumulators();
  pairs();
  if (failures) { printf("%d check(s) FAILED\n", failures); return 1; }
  printf("all fp128 checks passed\n");
  return 0;
}
