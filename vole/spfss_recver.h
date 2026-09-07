#ifndef SPFSS_RECVER_FP_H__
#define SPFSS_RECVER_FP_H__
#include <iostream>
#include "emp-tool/emp-tool.h"
#include "vole/twokeyprp.h"
#include "vole/fields/field_config.h"

using namespace emp;

template <typename IO, typename FP> 
class SpfssRecverFp {
public:
  block *ggm_tree = nullptr, *m = nullptr;
  block *ggm_tag = nullptr;  // GGM leaf tags (ring GGM-tree check only)
  bool *b = nullptr;
  std::size_t choice_pos, depth, leave_n;
  IO *io;
  FP share;

  SpfssRecverFp(IO *io, std::size_t depth_in) {
    this->io = io;
    this->depth = depth_in;
    this->leave_n = 1 << (depth_in - 1);
    m = new block[depth - 1];
    b = new bool[depth - 1];
  }

  ~SpfssRecverFp() {
    if(m != nullptr) delete[] m;
    if(b != nullptr) delete[] b;
    if(ggm_tree != nullptr) delete[] ggm_tree;
    if(ggm_tag != nullptr) delete[] ggm_tag;
  }

  std::size_t get_index() {
    choice_pos = 0;
    for (std::size_t i = 0; i < depth - 1; ++i) {
      choice_pos <<= 1;
      if (!b[i])
        choice_pos += 1;
    }
    return choice_pos;
  }

  // receive the message and reconstruct the tree
  // j: position of the secret, begins from 0
  template <typename OT> void recv(OT *ot, IO *io2, std::size_t s) {
    ot->recv(m, b, depth - 1, io2, s);
    share.recv(io2);
  }

  // receive the message and reconstruct the tree
  // j: position of the secret, begins from 0
  // delta2 only use low 64 bits
  // Reconstruct the punctured tree and write the reduced leaves (low/MAC lane;
  // value lane 0) into ggm_tree_mem. When FPS is block-sized (FP61x2) the tree
  // is reconstructed IN PLACE in the caller's buffer; otherwise a scratch tree
  // is used. The leaf loop fuses the conversion with a lazily reduced sum.
  template<typename FPS>
  void compute(FPS *ggm_tree_mem, FPS delta2_mac) {
    block *tree;
    if (sizeof(FPS) == sizeof(block)) {
      tree = reinterpret_cast<block *>(ggm_tree_mem);
    } else {
      ggm_tree = new block[leave_n];
      tree = ggm_tree;
    }
    tree_ = tree;
    ggm_tree_reconstruction(b, m);
    tree[choice_pos] = zero_block;

    FP nodes_sum;
    unsigned pending = 0;
    for (std::size_t i = 0; i < leave_n; ++i) {
      block leaf = tree[i];
      ggm_tree_mem[i].set_low_from_block(leaf);
      nodes_sum.add_raw(ggm_tree_mem[i].getLow());
      if (++pending == FP::lazy_adds) { nodes_sum.reduce(); pending = 0; }
    }
    if (pending) nodes_sum.reduce();
    nodes_sum = share + nodes_sum;
    nodes_sum = nodes_sum.negate();
    ggm_tree_mem[choice_pos] = delta2_mac + nodes_sum;
  }

  // Leaves of the most recent compute() (raw blocks for the ring tag check).
  block *tree_ = nullptr;

  void ggm_tree_reconstruction(bool *b, block *m) {
    block *ggm_tree = tree_;
    std::size_t to_fill_idx = 0;
    TwoKeyPRP prp(zero_block, makeBlock(0, 1));
    for (std::size_t i = 1; i < depth; ++i) {
      to_fill_idx = to_fill_idx * 2;
      ggm_tree[to_fill_idx] = ggm_tree[to_fill_idx + 1] = zero_block;
      if (b[i - 1] == false) {
        layer_recover(i, 0, to_fill_idx, m[i - 1], &prp);
        to_fill_idx += 1;
      } else
        layer_recover(i, 1, to_fill_idx + 1, m[i - 1], &prp);
    }
  }

  void layer_recover(std::size_t depth, std::size_t lr, std::size_t to_fill_idx, block sum,
                     TwoKeyPRP *prp) {
    block *ggm_tree = tree_;
    std::size_t layer_start = 0;
    std::size_t item_n = 1 << depth;
    block nodes_sum = zero_block;
    std::size_t lr_start = lr == 0 ? layer_start : (layer_start + 1);

    for (std::size_t i = lr_start; i < item_n; i += 2)
      nodes_sum = nodes_sum ^ ggm_tree[i];
    ggm_tree[to_fill_idx] = nodes_sum ^ sum;
    if (depth == this->depth - 1)
      return;
    for (int64_t i = (int64_t)item_n - 2; i >= 0; i -= 2)
      prp->node_expand_2to4(&ggm_tree[i * 2], &ggm_tree[i]);
  }

  // GGM-tree tag check (ring only; MozZ2karella Fig. 5 step 6). This party (P_S)
  // is punctured at choice_pos. Split into compute-only phases so MpfssReg can
  // batch all trees' messages together (see SpfssSenderFp).
  //
  // Phase B: tag every reconstructed leaf and recover the missing tag
  // t_{choice_pos} = K_top XOR (XOR_{j!=choice} t_j) from the sender's commitment.
  void ggm_tag_top(block K_top) {
    ggm_tag = new block[leave_n];
    ggm_leaf_tag(ggm_tag, tree_, leave_n);  // ggm_tag[choice_pos] is garbage
    block rest = zero_block;
    for (std::size_t i = 0; i < leave_n; ++i)
      if (i != choice_pos) rest = rest ^ ggm_tag[i];
    ggm_tag[choice_pos] = K_top ^ rest;
  }

  // Phase D: recompute Gamma' = <xi, t> over GF(2^128) for comparison against the
  // sender's Gamma (xi shared across all trees; no I/O).
  block ggm_tag_response(const block *xi) {
    block Gamma_prime;
    emp::vector_inn_prdt_sum_red(&Gamma_prime, ggm_tag, xi, (int)leave_n);
    return Gamma_prime;
  }

  // beta only use high 64 bits
  // x is the high 64 bits of z
  template<typename FPS>
  void consistency_check_msg_gen(FP &chi_alpha, FP &W,
                                 IO *io2, FPS beta, block seed,
                                 const FPS *input, std::size_t tree_idx) {
    if (field_traits<FP>::is_ring) {
      block tseed = ring_tree_seed(seed, tree_idx);
      bool *chi = new bool[leave_n];
      gen_binary_chi(chi, tseed, (int)leave_n);
      chi_alpha = chi[choice_pos] ? FP(1) : FP(0);
      FP v;
      v.setZero();
      for (std::size_t i = 0; i < leave_n; ++i)
        if (chi[i]) v = v + input[i].getLow();
      W = v;
      delete[] chi;
    } else {
      // Per-tree coefficients (paper Fig. 11: chi_i[j] independent across
      // trees i): derive a tree seed AES_seed(tree_idx) first, then the
      // polynomial-hash generator d = H(tseed), chi_j = d^(j+1). Both parties
      // hold the same (seed, tree_idx), so they agree.
      block tseed = ring_tree_seed(seed, tree_idx);
      FP *chi = new FP[leave_n];
      FP digest;
      digest.from_block(tseed);  // tseed = AES_seed(tree_idx) is already pseudorandom; no hash needed
      field_uni_hash_coeff_gen(chi, digest, leave_n);

      chi_alpha = chi[choice_pos];

      // W = \sum{chi_i*w_i}
      FP v;
      for(std::size_t i = 0; i < leave_n; ++i) {
        v = v + (chi[i] * input[i].getLow());
      }
      W = v;

      delete[] chi;
    }
  }
};
#endif
