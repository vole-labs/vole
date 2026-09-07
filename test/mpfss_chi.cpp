// Local check of the MPFSS consistency-check coefficients (paper Fig. 11):
//  * independence: with the same worker seed and identical leaf inputs, two
//    different trees must produce different check values V (chi_i differs
//    per tree i);
//  * agreement: for a receiver whose leaves equal the sender's, the receiver's
//    W equals the sender's V for every tree, and chi_alpha matches the
//    coefficient the sender used at the receiver's point.
// Runs the field path (fp61, fp107, f2k) and the ring path (z2k).
#include "emp-tool/emp-tool.h"
#include "vole/spfss_sender.h"
#include "vole/spfss_recver.h"
#include "vole/fields/fp61.h"
#include "vole/fields/fp61x2.h"
#include "vole/fields/fp107.h"
#include "vole/fields/fp107x2.h"
#include "vole/fields/fp2x128.h"
#include "vole/fields/z2k.h"

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
  const std::size_t depth = 6, leave_n = 1 << (depth - 1), trees = 8;
  PRG prg;
  block seed;
  prg.random_block(&seed, 1);

  // Common leaf values (as the sender's FP vector and the receiver's FPS
  // vector with the same low/MAC lane).
  std::vector<FP> leaves(leave_n);
  std::vector<FPS> leaves_vm(leave_n);
  for (std::size_t j = 0; j < leave_n; ++j) {
    leaves[j].rand(prg);
    leaves_vm[j].setLowHigh(leaves[j], FP(0));
  }

  std::vector<FP> V(trees);
  std::size_t distinct_pairs = 0, pairs = 0;
  for (std::size_t i = 0; i < trees; ++i) {
    SpfssSenderFp<NetIO, FP> snd(nullptr, depth, makeBlock(11, (uint64_t)i));
    V[i] = snd.consistency_check_msg_gen(nullptr, seed, leaves.data(), i);
  }
  for (std::size_t a = 0; a < trees; ++a)
    for (std::size_t b = a + 1; b < trees; ++b) {
      ++pairs;
      if (!(V[a] == V[b])) ++distinct_pairs;
    }
  printf("  distinct V across tree pairs: %zu / %zu\n", distinct_pairs, pairs);
  CHECK(distinct_pairs == pairs, "chi differs between trees (same seed, same leaves)");

  bool agree = true, alpha_ok = true;
  for (std::size_t i = 0; i < trees; ++i) {
    SpfssRecverFp<NetIO, FP> rcv(nullptr, depth);
    for (std::size_t l = 0; l < depth - 1; ++l) { bool b; prg.random_bool(&b, 1); rcv.b[l] = b; }
    rcv.get_index();
    FP chi_alpha, W;
    rcv.consistency_check_msg_gen(chi_alpha, W, nullptr, leaves_vm[0], seed,
                                  leaves_vm.data(), i);
    if (!(W == V[i])) agree = false;
    // Recompute the sender-side coefficient at choice_pos to compare.
    FP expect;
    if (field_traits<FP>::is_ring) {
      bool *chi = new bool[leave_n];
      gen_binary_chi(chi, ring_tree_seed(seed, i), (int)leave_n);
      expect = chi[rcv.choice_pos] ? FP(1) : FP(0);
      delete[] chi;
    } else {
      block tseed = ring_tree_seed(seed, i);
      FP d; d.from_block(tseed);  // same derivation as spfss_*: d = AES_seed(tree_idx)
      std::vector<FP> chi(leave_n);
      field_uni_hash_coeff_gen(chi.data(), d, (int)leave_n);
      expect = chi[rcv.choice_pos];
    }
    if (!(chi_alpha == expect)) alpha_ok = false;
  }
  CHECK(agree, "receiver W equals sender V for every tree");
  CHECK(alpha_ok, "receiver chi_alpha equals the sender's coefficient at alpha");
}

int main() {
  run<FP61, FP61x2>("fp61");
  run<FP107, FP107x2>("fp107");
  run<FP2x128, FP2x128x2>("f2k");
  run<Z2k64, Z2k64x2>("z2k");
  if (failures) { printf("%d check(s) FAILED\n", failures); return 1; }
  printf("all MPFSS chi checks passed\n");
  return 0;
}
