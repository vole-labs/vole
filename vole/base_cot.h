#ifndef VOLE_PROG_COT_H__
#define VOLE_PROG_COT_H__

// Correlated OTs for the GGM tree levels, on an emp-ot 1.0 OT extension
// (malicious-secure; default base OT CSW). Replaces the 0.3.0 BaseCot + IKNP
// pair. This is the only OT extension in the repo (COPE uses CSW base OTs
// directly). Engine: SoftSpoken<VOLE_COT_SOFTSPOKEN_K> (default k = 4), or
// IKNP with -DVOLE_COT_IKNP. Measured per 55k-COT round on one AWS box
// (localhost, malicious): IKNP 1.33 ms / 0.88 MB, SoftSpoken<2> 0.80 ms /
// 0.44 MB, SoftSpoken<4> 0.66 ms / 0.22 MB, SoftSpoken<8> 1.52 ms / 0.11 MB.
// The step is ~1% of an fp61 extend round's time and ~11% of its traffic
// (2.0 MB per round, mostly the GGM tree messages); k = 4 wins on both
// compute and traffic.
//
// Conventions kept from 0.3.0 (the F_2 instantiation relies on them):
//   * ot_delta has LSB 1 (main's OTExtension pins it);
//   * sender keys are stored with LSB 0, receiver blocks with LSB = choice;
//     (K xor b*Delta) & ~1 = (K & ~1) xor b*(Delta & ~1), so clearing on both
//     sides and re-inserting the bit keeps R = K' xor b*Delta exact.
// main's send_cot / recv_cot are random COT + 1-bit derandomisation, so the
// receiver's chosen bits are honoured but its raw blocks carry a random LSB;
// the wrapper rewrites the LSBs to the chosen bits.

#include "emp-ot/emp-ot.h"
#include "vole/preot.h"

template<typename IO>
class ProgBaseCot { public:
#ifdef VOLE_COT_IKNP
  using Engine = emp::IKNP;
#else
#ifndef VOLE_COT_SOFTSPOKEN_K
#define VOLE_COT_SOFTSPOKEN_K 4
#endif
  using Engine = emp::SoftSpoken<VOLE_COT_SOFTSPOKEN_K>;
#endif
  int party;
  IO *io;
  bool malicious;
  emp::block one, minusone;
  emp::block ot_delta;
  Engine *ext = nullptr;

  ProgBaseCot(int party, IO *io, bool malicious = true)
      : party(party), io(io), malicious(malicious) {
    minusone = emp::makeBlock(0xFFFFFFFFFFFFFFFFLL, 0xFFFFFFFFFFFFFFFELL);
    one = emp::makeBlock(0x0LL, 0x1LL);
    ext = new Engine(party, io, malicious);
  }
  ~ProgBaseCot() { delete ext; }

  // Base OTs (0.3.0's cot_gen_pre): main bootstraps lazily on first use.
  void cot_gen_pre() {
    if (party == emp::ALICE) ot_delta = ext->Delta;
  }

  // ALICE: keys for `size` COTs into pre_ot (LSB cleared).
  void prog_cot_gen(emp::OTPre<IO> *pre_ot, std::size_t size) {
    emp::block *ot_data = new emp::block[size];
    ext->send_cot(ot_data, (int64_t)size);
    io->flush();
    ot_delta = ext->Delta;
    for (std::size_t i = 0; i < size; ++i) ot_data[i] = ot_data[i] & minusone;
    pre_ot->send_pre(ot_data, ot_delta);
    delete[] ot_data;
  }

  // BOB: blocks for the chosen bits pre_bool_ini (LSB rewritten to the bit).
  void prog_cot_gen(emp::OTPre<IO> *pre_ot, bool *pre_bool_ini, std::size_t size) {
    emp::block *ot_data = new emp::block[size];
    ext->recv_cot(ot_data, pre_bool_ini, (int64_t)size);
    io->flush();
    emp::block ch[2] = {emp::zero_block, one};
    for (std::size_t i = 0; i < size; ++i)
      ot_data[i] = (ot_data[i] & minusone) ^ ch[pre_bool_ini[i]];
    pre_ot->recv_pre(ot_data, pre_bool_ini);
    delete[] ot_data;
  }
};

#endif  // VOLE_PROG_COT_H__
