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

  // send the nodes by oblivious transfer
  void compute(FP *ggm_tree_mem, FP secret,
               FP gamma_mac) {
    this->delta = secret;
    ggm_tree = new block[leave_n];
    ggm_tree_gen(m, m + depth - 1);

    secret_sum.setZero();
    for (std::size_t i = 0; i < leave_n; ++i) {
      ggm_tree_mem[i].from_block(ggm_tree[i]);
      secret_sum = secret_sum + ggm_tree_mem[i];
    }
    FP zz(0);
    if (FP::PR_num_pack > 1) {
      for (std::size_t i = 0; i < leave_n; ++i) {
        ggm_tree_mem[i].setHigh(zz);
      }
    }
    secret_sum = secret_sum.negate();
    secret_sum = gamma_mac + secret_sum;
  }

  // send the nodes by oblivious transfer
  template <typename OT> void send(OT *ot, IO *io2, std::size_t s) {
    ot->send(m, &m[depth - 1], depth - 1, io2, s);
    secret_sum.send(io2);
    io2->flush();
  }

  // generate GGM tree from the top
  void ggm_tree_gen(block *ot_msg_0, block *ot_msg_1) {
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
    ggm_leaf_tag(ggm_tag, ggm_tree, leave_n);
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
