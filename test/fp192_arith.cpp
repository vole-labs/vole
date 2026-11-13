// Local (no network) checks of the p192 field, p = 2^192 - 2^64 - 1
// (vole/fields/fp192.h, fp192x2.h):
//  1) known answers: +, -, *, inv against Python big-int results (fp192_kat.h),
//     including the carry/borrow edge values;
//  2) reduction edges: the non-canonical inputs [p, 2^192) are reduced by the
//     U192 constructor, rand and from_2blocks; 2^192 = 2^64 + 1 (mod p);
//  3) algebraic identities on random elements (associativity, distributivity,
//     inverse, negation, shifts, the int / uint64 constructors);
//  4) the lazy accumulators FP192::Acc / FP192x2::Acc used by the LPN encoders;
//  5) FP192x2 lane-wise operations agree with FP192.
//
// usage: test_fp192_arith
#include "emp-tool/emp-tool.h"
#include "vole/fields/fp192.h"
#include "vole/fields/fp192x2.h"
#include "test/fp192_kat.h"
#include <string>

using namespace emp;

static int failures = 0;
#define CHECK(cond, msg)                                                   \
  do {                                                                     \
    if (!(cond)) { printf("  FAIL: %s\n", msg); ++failures; }              \
    else         { printf("  ok  : %s\n", msg); }                          \
  } while (0)

static U192 parse(const char *hex) {  // "0x" + 48 hex digits
  std::string s(hex + 2);
  U192 r{{0, 0, 0}};
  for (int limb = 0; limb < 3; ++limb)
    r.w[limb] = std::stoull(s.substr(48 - 16 * (limb + 1), 16), nullptr, 16);
  return r;
}
static FP192 raw(const U192 &v) { FP192 r; r.val = v; return r; }  // no reduction

int main() {
  printf("[p192 known answers]\n");
  {
    std::size_t n = sizeof(KAT) / sizeof(KAT[0]), bad_add = 0, bad_sub = 0,
                bad_mul = 0, bad_inv = 0, bad_neg = 0;
    for (std::size_t i = 0; i < n; ++i) {
      FP192 a = raw(parse(KAT[i][0])), b = raw(parse(KAT[i][1]));
      if (a + b != raw(parse(KAT[i][2]))) ++bad_add;
      if (a - b != raw(parse(KAT[i][3]))) ++bad_sub;
      if (a * b != raw(parse(KAT[i][4]))) ++bad_mul;
      if (a.inv() != raw(parse(KAT[i][5]))) ++bad_inv;
      if (b - a != -(a - b)) ++bad_neg;
    }
    printf("  %zu vectors\n", n);
    CHECK(bad_add == 0, "a + b");
    CHECK(bad_sub == 0, "a - b");
    CHECK(bad_mul == 0, "a * b");
    CHECK(bad_inv == 0, "inv(a) (0 -> 0)");
    CHECK(bad_neg == 0, "b - a == -(a - b)");
  }

  printf("[p192 reduction edges]\n");
  {
    const U192 P = FP192::PR;
    CHECK(FP192(P) == FP192(0), "U192(p) reduces to 0");
    U192 pm1 = P; pm1.w[0] -= 1;
    CHECK(FP192(pm1).val == pm1, "U192(p-1) is kept");
    U192 top{{~0ULL, ~0ULL, ~0ULL}};  // 2^192 - 1 = p + 2^64
    CHECK(FP192(top) == FP192(U192{{0, 1, 0}}), "U192(2^192 - 1) reduces to 2^64");
    // 2^192 = 2^191 * 2 = 2^64 + 1
    FP192 two191 = FP192(U192{{0, 0, 1ULL << 63}});
    CHECK(two191 + two191 == FP192(U192{{1, 1, 0}}), "2^191 + 2^191 = 2^64 + 1");
    CHECK((FP192(1) << 191) * FP192(2) == FP192(U192{{1, 1, 0}}), "(1 << 191) * 2 = 2^64 + 1");
    // from_2blocks: 192 bits taken as (b0, low half of b1), reduced.
    block ones = makeBlock(~0ULL, ~0ULL);
    FP192 f; f.from_2blocks(ones, ones);
    CHECK(f == FP192(U192{{0, 1, 0}}), "from_2blocks(all ones) = 2^64");
    block lo = makeBlock(0x1111111111111111ULL, 0x2222222222222222ULL);
    block hi = makeBlock(0xdeadbeefULL, 0x3333333333333333ULL);
    f.from_2blocks(lo, hi);
    CHECK(f.val == (U192{{0x2222222222222222ULL, 0x1111111111111111ULL, 0x3333333333333333ULL}}),
          "from_2blocks limb order (b0 low, b0 high, b1 low)");
    FP192 g; g.from_block(lo);
    CHECK(g.val == (U192{{0x2222222222222222ULL, 0x1111111111111111ULL, 0}}),
          "from_block is the block's integer value");
    // rand / recv-style canonicalisation: every sample < p.
    PRG prg;
    bool canon = true;
    for (int i = 0; i < 100000; ++i) { FP192 r; r.rand(prg); if (FP192::geq_p(r.val)) canon = false; }
    CHECK(canon, "rand() always returns a canonical value");
    CHECK(FP192(-1) == raw(pm1), "FP192(-1) = p - 1");
    CHECK(FP192(-5) + FP192(5) == FP192(0), "FP192(-5) + 5 = 0");
    CHECK(FP192(0).inv() == FP192(0), "inv(0) = 0");
  }

  printf("[p192 identities on random elements]\n");
  {
    PRG prg;
    bool assoc = true, distr = true, inv = true, neg = true, shift = true, u64 = true;
    for (int i = 0; i < 20000; ++i) {
      FP192 a, b, c;
      a.rand(prg); b.rand(prg); c.rand(prg);
      if ((a * b) * c != a * (b * c)) assoc = false;
      if (a * (b + c) != a * b + a * c) distr = false;
      if (i < 2000 && !(a == FP192(0)) && a * a.inv() != FP192(1)) inv = false;
      if (a + (-a) != FP192(0) || a.negate() != -a) neg = false;
      int n = i % 192;
      FP192 pw(1);
      for (int k = 0; k < n; ++k) pw = pw + pw;
      if ((a << n) != a * pw) shift = false;
      uint64_t s; prg.random_data(&s, sizeof(s));
      if (a * s != a * FP192(s) || a + s != a + FP192(s)) u64 = false;
    }
    CHECK(assoc, "(ab)c = a(bc)");
    CHECK(distr, "a(b + c) = ab + ac");
    CHECK(inv, "a * inv(a) = 1");
    CHECK(neg, "a + (-a) = 0, negate() = -a");
    CHECK(shift, "a << n = a * 2^n for n < 192");
    CHECK(u64, "uint64 operand overloads");
  }

  printf("[p192 lazy accumulator (LPN encoders)]\n");
  {
    PRG prg;
    bool ok = true, ok2 = true;
    U192 pm1 = FP192::PR; pm1.w[0] -= 1;
    for (int n : {1, 2, 10, 11, 100, 5000}) {
      for (int rep = 0; rep < 20; ++rep) {
        FP192 x0; x0.rand(prg);
        if (rep == 0) x0 = raw(pm1);
        FP192::Acc acc(x0);
        FP192 ref = x0;
        FP192x2 y0(x0, x0 + FP192(1));
        FP192x2::Acc acc2(y0);
        FP192x2 ref2 = y0;
        for (int i = 0; i < n; ++i) {
          FP192 y; y.rand(prg);
          if (rep == 0) y = raw(pm1);  // worst case: every add carries
          acc.add(y); ref = ref + y;
          FP192x2 z(y, -y); acc2.add(z); ref2 = ref2 + z;
        }
        if (acc.get() != ref) ok = false;
        if (acc2.get() != ref2) ok2 = false;
      }
    }
    CHECK(ok, "FP192::Acc equals the reduced running sum (1..5000 terms, p-1 worst case)");
    CHECK(ok2, "FP192x2::Acc equals the reduced running sum");
  }

  printf("[p192 FP192x2 lanes]\n");
  {
    PRG prg;
    bool ok = true;
    for (int i = 0; i < 1000; ++i) {
      FP192 a, b, c, d, s;
      a.rand(prg); b.rand(prg); c.rand(prg); d.rand(prg); s.rand(prg);
      FP192x2 x(a, b), y(c, d);
      if ((x + y) != FP192x2(a + c, b + d)) ok = false;
      if ((x - y) != FP192x2(a - c, b - d)) ok = false;
      if ((x * y) != FP192x2(a * c, b * d)) ok = false;
      if ((x * s) != FP192x2(a * s, b * s)) ok = false;
      if ((x + s) != FP192x2(a + s, b + s)) ok = false;
      if (x.negate() != FP192x2(-a, -b)) ok = false;
      FP192x2 z = x; z.add_raw(y); z.reduce();
      if (z != x + y) ok = false;
      if (x.getLow() != a || x.getHigh() != b) ok = false;
    }
    CHECK(ok, "lane-wise +, -, *, scalar ops, negate, add_raw, getLow/getHigh");
  }

  if (failures) { printf("%d check(s) FAILED\n", failures); return 1; }
  printf("all p192 field checks passed\n");
  return 0;
}
