#ifndef SPFSS_SENDER_FP_H__
#define SPFSS_SENDER_FP_H__
#include <iostream>
#include "emp-ot/emp-ot.h"
#include "emp-tool/emp-tool.h"
#include "vole/fields/field_config.h"

using namespace emp;

template <typename IO, typename FP> 
class SpfssSenderFp {
public:
  block seed;
  block *ggm_tree = nullptr, *m = nullptr;
  block *ggm_tag = nullptr;  // GGM leaf tags (ring GGM-tree check only)
  FP delta;
  FP secret_sum;
  IO *io;
  std::size_t depth;
  std::size_t leave_n;
  PRG prg;

  SpfssSenderFp(IO *io, std::size_t depth_in) {
    initialization(io, depth_in);
    prg.random_block(&seed, 1);
  }

  void initialization(IO *io, std::size_t depth_in) {
    this->io = io;
    this->depth = depth_in;
    this->leave_n = 1 << (this->depth - 1);
    m = new block[(depth - 1) * 2];
  }

  ~SpfssSenderFp() {
    if(m != nullptr) delete[] m;
    if(ggm_tree != nullptr) delete[] ggm_tree;
    if(ggm_tag != nullptr) delete[] ggm_tag;
  }

  // Expand the GGM tree and write the reduced leaves into ggm_tree_mem. When
  // FP is block-sized the tree is expanded IN PLACE in the caller's buffer
  // (no scratch allocation, no copy); otherwise (8-byte FP61) a scratch tree
  // is used. The leaf loop fuses the field conversion with a lazily reduced
  // sum.
  void compute(FP *ggm_tree_mem, FP secret,
               FP gamma_mac) {
    this->delta = secret;
    block *tree;
    if (sizeof(FP) == sizeof(block)) {
      tree = reinterpret_cast<block *>(ggm_tree_mem);
    } else {
      ggm_tree = new block[leave_n];
      tree = ggm_tree;
    }
    tree_ = tree;
    ggm_tree_gen(m, m + depth - 1);

    secret_sum.setZero();
    unsigned pending = 0;
    for (std::size_t i = 0; i < leave_n; ++i) {
      block leaf = tree[i];
      ggm_tree_mem[i].from_block(leaf);
      secret_sum.add_raw(ggm_tree_mem[i]);
      if (++pending == FP::lazy_adds) { secret_sum.reduce(); pending = 0; }
    }
    if (pending) secret_sum.reduce();
    secret_sum = secret_sum.negate();
    secret_sum = gamma_mac + secret_sum;
  }

  // Leaves of the most recent compute() (raw blocks for the ring tag check;
  // identical to the field values when from_block is the identity).
  block *tree_ = nullptr;

  // send the nodes by oblivious transfer
  template <typename OT> void send(OT *ot, IO *io2, std::size_t s) {
    ot->send(m, &m[depth - 1], depth - 1, io2, s);
    secret_sum.send(io2);
    io2->flush();
  }

  // generate GGM tree from the top
  void ggm_tree_gen(block *ot_msg_0, block *ot_msg_1) {
    block *ggm_tree = tree_;
    TwoKeyPRP *prp = new TwoKeyPRP(zero_block, makeBlock(0, 1));
    prp->node_expand_1to2(ggm_tree, seed);
    ot_msg_0[0] = ggm_tree[0];
    ot_msg_1[0] = ggm_tree[1];
    for (std::size_t h = 1; h < depth - 1; ++h) {
      ot_msg_0[h] = ot_msg_1[h] = zero_block;
      std::size_t sz = 1 << h;
      for (int64_t i = (int64_t)sz - 2; i >= 0; i -= 2) {
        prp->node_expand_2to4(&ggm_tree[i * 2], &ggm_tree[i]);
        ot_msg_0[h] = ot_msg_0[h] ^ ggm_tree[i * 2];
        ot_msg_0[h] = ot_msg_0[h] ^ ggm_tree[i * 2 + 2];
        ot_msg_1[h] = ot_msg_1[h] ^ ggm_tree[i * 2 + 1];
        ot_msg_1[h] = ot_msg_1[h] ^ ggm_tree[i * 2 + 3];
      }
    }
    delete prp;
  }

  // GGM-tree tag check (ring only; MozZ2karella Fig. 5 step 6). This party (P_R)
  // holds the full tree. Split into compute-only phases so MpfssReg can batch the
  // three messages across ALL trees (one K_top buffer, one shared challenge, one
  // Gamma buffer) instead of a per-tree ping-pong.
  //
  // Phase A: tag every leaf and commit K_top = XOR_j t_j (no I/O).
  block ggm_tag_top() {
    ggm_tag = new block[leave_n];
    ggm_leaf_tag(ggm_tag, tree_, leave_n);
    block K_top = zero_block;
    for (std::size_t i = 0; i < leave_n; ++i)
      K_top = K_top ^ ggm_tag[i];
    return K_top;
  }

  // Phase C: answer the universal-hash challenge with Gamma = <xi, t> over
  // GF(2^128) (xi is shared across all trees in the MPFSS; no I/O).
  block ggm_tag_response(const block *xi) {
    block Gamma;
    emp::vector_inn_prdt_sum_red(&Gamma, ggm_tag, xi, (int)leave_n);
    return Gamma;
  }

  // consistency check: Protocol PI_spsVOLE.
  // Field path: V = sum_i chi_i*v_i with random ring/field coefficients chi.
  // Ring path (Z_{2^k}): binary chi in {0,1}^n of weight n/2 (MozZ2karella Fig. 5
  // step 8), per-tree-seeded so coefficients differ across trees.
  FP consistency_check_msg_gen(IO *io2, block seed, const FP *input,
                               std::size_t tree_idx) {
    FP V;
    if (field_traits<FP>::is_ring) {
      block tseed = ring_tree_seed(seed, tree_idx);
      bool *chi = new bool[leave_n];
      gen_binary_chi(chi, tseed, (int)leave_n);
      V.setZero();
      for (std::size_t i = 0; i < leave_n; ++i)
        if (chi[i]) V = V + input[i];
      delete[] chi;
    } else {
      // Per-tree coefficients (paper Fig. 11: chi_i[j] independent across
      // trees i): derive a tree seed AES_seed(tree_idx) first, then the
      // polynomial-hash generator d = H(tseed), chi_j = d^(j+1). Both parties
      // hold the same (seed, tree_idx), so they agree.
      block tseed = ring_tree_seed(seed, tree_idx);
      FP *chi = new FP[leave_n];
      FP digest;
      digest.from_block(Hash::hash_for_block(&tseed, sizeof(block)));
      uni_hash_coeff_gen(chi, digest, leave_n);
      // V = \sum{chi_i*v_i}
      V = vector_inn_prdt_sum_red(chi, input, leave_n);
      delete[] chi;
    }
    return V;
  }
};

#endif
