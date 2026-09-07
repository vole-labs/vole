#ifndef VOLE_AES_NI_H__
#define VOLE_AES_NI_H__

// Plain AES-NI kernels for this repo's two hot loops (GGM node expansion and
// the LPN index PRP). emp-tool 1.0 routes AES through a tiered VAES/AES-NI
// dispatch (EMP_AES_TARGET_ATTR functions) that is tuned for large batches;
// for the 2- to 10-block calls made millions of times per extend round it
// stopped inlining into our loops and cost 1.6x (LPN) to 3x (GGM). These are
// the 0.3.0 kernels (10-round AES-128 on emp's expanded AES_KEY), fully
// inline so the compiler pipelines them across consecutive calls.
//
// Only the key schedule (AES_set_encrypt_key / AES_KEY::rd_key[11]) is taken
// from emp-tool; large batches (ComMatrixFp, CCRH) keep using emp-tool.

#include "emp-tool/emp-tool.h"

namespace vole {

// blks[i] <- AES_key(blks[i]) for i < nblks, all rounds interleaved.
inline void aes_ecb_encrypt_blks(emp::block *blks, int nblks, const emp::AES_KEY *key) {
  for (int i = 0; i < nblks; ++i) blks[i] = _mm_xor_si128(blks[i], key->rd_key[0]);
  for (int j = 1; j < 10; ++j)
    for (int i = 0; i < nblks; ++i) blks[i] = _mm_aesenc_si128(blks[i], key->rd_key[j]);
  for (int i = 0; i < nblks; ++i) blks[i] = _mm_aesenclast_si128(blks[i], key->rd_key[10]);
}

// numKeys keys x numEncs blocks each, laid out blks[k * numEncs + e]
// (emp's ParaEnc<numKeys, numEncs> layout), all interleaved.
template <int numKeys, int numEncs>
inline void para_enc(emp::block *blks, const emp::AES_KEY *keys) {
  for (int k = 0; k < numKeys; ++k)
    for (int e = 0; e < numEncs; ++e)
      blks[k * numEncs + e] = _mm_xor_si128(blks[k * numEncs + e], keys[k].rd_key[0]);
  for (int j = 1; j < 10; ++j)
    for (int k = 0; k < numKeys; ++k)
      for (int e = 0; e < numEncs; ++e)
        blks[k * numEncs + e] = _mm_aesenc_si128(blks[k * numEncs + e], keys[k].rd_key[j]);
  for (int k = 0; k < numKeys; ++k)
    for (int e = 0; e < numEncs; ++e)
      blks[k * numEncs + e] = _mm_aesenclast_si128(blks[k * numEncs + e], keys[k].rd_key[10]);
}

}  // namespace vole

#endif  // VOLE_AES_NI_H__
