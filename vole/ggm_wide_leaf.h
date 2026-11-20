#ifndef VOLE_GGM_WIDE_LEAF_H__
#define VOLE_GGM_WIDE_LEAF_H__
// Leaf expansion for fields wider than a block (wide_leaf<FP>), shared by the
// SPFSS sender and receiver so both derive identical elements.
//
//   2 blocks (FP192, FP25519): leaf -> (c0, c1), one more GGM level.
//   3 blocks (FPP256):         leaf -> (c0, c1), c1 -> (d0, d1); use
//                              (c0, d0, d1): two more levels down the right
//                              child, a PRG of the leaf like the rest of the
//                              tree.
// Leaves are processed four at a time (node_expand_4to8), then one at a time.
// emit(i, b) receives leaf i's leaf_blocks<FP> blocks.

#include "emp-tool/emp-tool.h"
#include "vole/twokeyprp.h"

template <int NB, typename Emit>
void expand_wide_leaves(const block *tree, std::size_t leave_n, Emit emit) {
  static_assert(NB == 2 || NB == 3, "wide leaves are two or three blocks");
  TwoKeyPRP prp(zero_block, makeBlock(0, 1));
  block e1[8], e2[8];
  std::size_t i = 0;
  for (; i + 4 <= leave_n; i += 4) {
    block par[4] = {tree[i], tree[i + 1], tree[i + 2], tree[i + 3]};
    prp.node_expand_4to8(e1, par);  // children of leaf i+k: e1[2k], e1[2k+1]
    if constexpr (NB == 2) {
      for (int k = 0; k < 4; ++k) {
        const block b[2] = {e1[2 * k], e1[2 * k + 1]};
        emit(i + k, b);
      }
    } else {
      block odd[4] = {e1[1], e1[3], e1[5], e1[7]};
      prp.node_expand_4to8(e2, odd);  // children of e1[2k+1]: e2[2k], e2[2k+1]
      for (int k = 0; k < 4; ++k) {
        const block b[3] = {e1[2 * k], e2[2 * k], e2[2 * k + 1]};
        emit(i + k, b);
      }
    }
  }
  for (; i < leave_n; ++i) {
    prp.node_expand_1to2(e1, tree[i]);
    if constexpr (NB == 2) {
      const block b[2] = {e1[0], e1[1]};
      emit(i, b);
    } else {
      prp.node_expand_1to2(e2, e1[1]);
      const block b[3] = {e1[0], e2[0], e2[1]};
      emit(i, b);
    }
  }
}

#endif  // VOLE_GGM_WIDE_LEAF_H__
