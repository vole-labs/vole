// Local (no network) check of the single-point FSS leaf conversion
// (SpfssSenderFp::compute / SpfssRecverFp::compute), every field:
//  * the receiver, given the OT messages for its choice bits, reproduces the
//    sender's leaf value at every position except alpha;
//  * at alpha, w[alpha] - v[alpha] = delta2 - gamma (the sender's gamma_mac and
//    the receiver's delta2_mac), i.e. the punctured slot is consistent with the
//    sender's published sum;
//  * for fields wider than a block (p192) every leaf is a
//    full-width element: with a single 128-bit block per leaf the top limb would
//    always be zero and v[alpha] = w[alpha] - beta*Delta would leak Delta.
// Tree depths 2 and 3 exercise the 1-leaf tail of the 4-leaf batched expansion
// (leave_n = 2, 4); depth 9 is a full tree of 256 leaves.
//
// usage: test_ggm_leaf
#include "emp-tool/emp-tool.h"
#include "vole/spfss_sender.h"
#include "vole/spfss_recver.h"
#include "vole/fields/field_config.h"

using namespace emp;

static int failures = 0;
#define CHECK(cond, msg)                                                   \
  do {                                                                     \
    if (!(cond)) { printf("  FAIL: %s\n", msg); ++failures; }              \
    else         { printf("  ok  : %s\n", msg); }                          \
  } while (0)

template <typename FP, typename FPS>
void run(const char *tag) {
  printf("[%s]\n", tag);
  PRG prg;
  bool others = true, alpha = true, wide = true;
  std::size_t trials = 0;
  for (std::size_t depth : {2, 3, 9}) {
    const std::size_t leave_n = 1 << (depth - 1);
    for (int rep = 0; rep < 8; ++rep, ++trials) {
      block seed;
      prg.random_block(&seed, 1);
      FP gamma, secret;
      gamma.rand(prg);
      secret.rand(prg);
      FPS delta2;
      { FP lo, hi; lo.rand(prg); hi.rand(prg); delta2.setLowHigh(lo, hi); }

      std::vector<FP> v(leave_n);
      SpfssSenderFp<NetIO, FP> snd(nullptr, depth, seed);
      snd.compute(v.data(), secret, gamma);

      SpfssRecverFp<NetIO, FP> rcv(nullptr, depth);
      for (std::size_t l = 0; l < depth - 1; ++l) {
        bool b; prg.random_bool(&b, 1);
        rcv.b[l] = b;
        rcv.m[l] = b ? snd.m[(depth - 1) + l] : snd.m[l];  // OT: choice b gets msg_b
      }
      rcv.get_index();
      rcv.share = snd.secret_sum;
      std::vector<FPS> w(leave_n);
      rcv.compute(w.data(), delta2);

      for (std::size_t i = 0; i < leave_n; ++i)
        if (i != rcv.choice_pos && w[i].getLow() != v[i]) others = false;
      if (w[rcv.choice_pos].getLow() - v[rcv.choice_pos] != delta2.getLow() - gamma)
        alpha = false;
      if constexpr (FP::PR_bit_len > 128) {  // independent of the wide_leaf trait
        for (std::size_t i = 0; i < leave_n; ++i)
          if (v[i].val.w[2] == 0) wide = false;  // probability 2^-64 per leaf
      }
    }
  }
  printf("  %zu trees\n", trials);
  CHECK(others, "receiver leaves equal the sender's off alpha");
  CHECK(alpha, "w[alpha] - v[alpha] = delta2 - gamma");
  if constexpr (FP::PR_bit_len > 128)
    CHECK(wide, "every leaf uses the full 192-bit width (two blocks per leaf)");
  CHECK(wide_leaf<FP>::value == (FP::PR_bit_len > 128 && !field_traits<FP>::is_ring),
        "wide_leaf set exactly for fields wider than a block");
}

int main() {
  run<FP61, FP61x2>("fp61");
  run<FP107, FP107x2>("fp107");
  run<FP192, FP192x2>("p192");
  run<FP2x128, FP2x128x2>("f2k");
  run<Z2k64, Z2k64x2>("z2k");
  if (failures) { printf("%d check(s) FAILED\n", failures); return 1; }
  printf("all GGM leaf checks passed\n");
  return 0;
}
