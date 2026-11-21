// Two-party check of the MPFSS final check (vole/mpfss_reg.h): the
// Delta-holder (ALICE) commits to vb first, the receiver (BOB) reveals va,
// and ALICE opens only if va == vb. One VoleTriple setup (whose stages run
// the check honestly) and one extend round, whose check runs in `mode`:
//
//   honest     both run the protocol                  -> both succeed
//   xstar      BOB sends x* + 1                        -> ALICE: "MPFSS batch check fails"
//                                                         (va != vb; no opening is sent)
//   noncanon   BOB sends x* + p (FP61 encoding >= p)   -> ALICE: "non-canonical value"
//   va         BOB reveals va + 1                      -> ALICE: "MPFSS batch check fails"
//   badopen    ALICE opens with another salt           -> BOB: "opening does not match"
//
// The expected party and message are checked by the CMake test (grep on its
// output); the extend round prints "round passed" only on success.
//
// usage: test_mpfss_check <party> <port> <mode> [field]   field in {fp61, p192, p256}
#include "emp-tool/emp-tool.h"
#include "vole/vole-fp.h"
#include "vole/fields/field_config.h"
#include "test.h"
#include <string>

using namespace emp;

int party, port;

template <typename FP, typename FPS>
void run(NetIO **ios, const std::string &mode, const char *tag) {
  VoleTriple<NetIO, FP, FPS> vt(party, 1, ios, ferret_b9);
  PRG prg;
  if (party == ALICE) {
    FP Delta;
    Delta.rand(prg);
    vt.template setup<FP>(Delta);
  } else {
    block s;
    prg.random_block(&s, 1);
    vt.setup_prog(s);
    vt.template setup<FP>();
  }
  printf("[%s] setup passed\n", tag);
  fflush(stdout);

  using M = MpfssRegFp<NetIO, FP, FPS>;
  if (party == BOB) {
    if (mode == "xstar") vt.mpfss->test_cheat = M::kCheatXStar;
    if (mode == "noncanon") vt.mpfss->test_cheat = M::kCheatNonCanonical;
    if (mode == "va") vt.mpfss->test_cheat = M::kCheatVa;
  } else if (mode == "badopen") {
    vt.mpfss->test_cheat = M::kCheatBadOpening;
  }

  if (party == ALICE) {
    std::vector<FP> buf(vt.param.n);
    vt.extend_send(buf.data());
  } else {
    std::vector<FPS> buf(vt.param.n);
    vt.extend_recv(buf.data());
  }
  printf("[%s] round passed\n", tag);
  fflush(stdout);
}

int main(int argc, char **argv) {
  parse_party_and_port(argv, &party, &port);
  std::string mode = argc > 3 ? argv[3] : "honest";
  std::string field = argc > 4 ? argv[4] : "fp61";
  if (mode != "honest" && mode != "xstar" && mode != "noncanon" && mode != "va" && mode != "badopen")
    error("unknown mode: honest | xstar | noncanon | va | badopen");
  NetIO *ios[1] = {new NetIO(party == ALICE ? nullptr : "127.0.0.1", port, true)};
  if (field == "fp61") run<FP61, FP61x2>(ios, mode, "fp61");
  else if (field == "p192") run<FP192, FP192x2>(ios, mode, "p192");
  else if (field == "p256") run<FPP256, FPP256x2>(ios, mode, "p256");
  else error("unknown field: fp61 | p192 | p256");
  delete ios[0];
  return 0;
}
