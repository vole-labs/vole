#include "emp-tool/emp-tool.h"
#include "vole/cvole.h"
// All supported fields compiled in; chosen at run time (argv[3]).
#include "vole/fields/fp61.h"
#include "vole/fields/fp61x2.h"
#include "vole/fields/fp107.h"
#include "vole/fields/fp107x2.h"
#include "vole/fields/fp2x128.h"
#include "vole/fields/z2k.h"
#include "test.h"

using namespace emp;
using namespace std;

int party, port;
const std::size_t threads = 1;

// Field-generic committed-VOLE test:
//  1) verify the VOLE relation M = K + value*Delta on x and on com,
//  2) run the committed-VOLE consistency check (publish com + Hash(M[com])).
template <typename FQ, typename FP, typename FPS>
void test_cvole(NetIO *ios[threads], int party, const char *tag) {
  // Small param + multiple rounds, to exercise self-seeding across rounds.
  const std::size_t rounds = 2;
  CVoleFp<NetIO, FP, FPS> cvole(party, threads, ios, cvole_fp_n2to14);

  FP Delta;
  PRG prg;

  block committer_seed = makeBlock(0xC0FFEEULL, 0x1234ULL);
  if (party == ALICE) {
    Delta.rand(prg);
    if (FP::PR_num_pack > 1)
      Delta.setHigh(FP(0));
    restrict_delta(Delta);  // ring: clamp Delta to its s-bit sub-ring (no-op for fields)
    cvole.setup(Delta);
  } else {
    cvole.setup_prog(committer_seed);  // before setup: reproducible bootstrap
    cvole.setup();
  }

  std::size_t xn = cvole.x_buf_size(rounds);    // overlapping x buffer
  std::size_t cn = cvole.com_buf_size(rounds);  // total |com|

  auto start = clock_start();
  if (party == ALICE) {
    FP *x_key = new FP[xn];
    FP *com_key = new FP[cn];
    std::size_t usable = cvole.extend_inplace_send(x_key, com_key, rounds);
    std::cout << "[" << tag << "] extend " << time_from(start) / 1000 << " ms, "
              << usable << " usable of " << xn << std::endl;
    if (usable != cvole.x_usable(rounds)) error("usable count mismatch");
    for (std::size_t i = usable; i < xn; ++i)   // reserved tail is scrubbed
      if (!(x_key[i] == FP(0))) error("reserved tail not zeroed (sender)");

    cvole.check_vole_send(x_key, xn);    // M[x]   = K[x]   + x*Delta
    cvole.check_vole_send(com_key, cn);  // M[com] = K[com] + com*Delta
    cvole.consistency_check_verify(com_key, cn);
    std::cout << "[" << tag << "] consistency check passed (rounds=" << rounds
              << ")" << std::endl;

    delete[] x_key;
    delete[] com_key;
  } else {
    FPS *x_vm = new FPS[xn];
    FPS *com_vm = new FPS[cn];
    std::size_t usable = cvole.extend_inplace_recv(x_vm, com_vm, rounds);
    std::cout << "[" << tag << "] extend " << time_from(start) / 1000 << " ms, "
              << usable << " usable of " << xn << std::endl;
    for (std::size_t i = usable; i < xn; ++i)   // reserved tail is scrubbed
      if (!(x_vm[i].getHigh() == FP(0)) || !(x_vm[i].getLow() == FP(0)))
        error("reserved tail not zeroed (receiver)");

    cvole.check_vole_recv(x_vm, xn);
    cvole.check_vole_recv(com_vm, cn);

    FP *com_val = new FP[cn];
    for (std::size_t i = 0; i < cn; ++i) com_val[i] = com_vm[i].getHigh();
    cvole.consistency_check_commit(com_val, com_vm, cn);
    std::cout << "[" << tag << "] consistency check passed (rounds=" << rounds
              << ")" << std::endl;

    // King-side local reproduction: a non-interactive instance with the same
    // seed must reproduce the committer's x and com VALUES exactly.
    CVoleFp<NetIO, FP, FPS> king(threads, cvole_fp_n2to14);
    king.setup_prog(committer_seed);
    king.setup_local();
    FP *x_loc = new FP[xn];
    FP *com_loc = new FP[cn];
    king.extend_inplace_local(x_loc, com_loc, rounds);
    for (std::size_t i = 0; i < xn; ++i)
      if (x_loc[i] != x_vm[i].getHigh()) error("king-local x mismatch");
    for (std::size_t i = 0; i < cn; ++i)
      if (com_loc[i] != com_vm[i].getHigh()) error("king-local com mismatch");
    if (com_loc[0] != com_val[0]) error("king-local com vs published mismatch");
    std::cout << "[" << tag << "] king-local x/com match committer" << std::endl;

    delete[] x_loc;
    delete[] com_loc;
    delete[] com_val;
    delete[] x_vm;
    delete[] com_vm;
  }
}

int main(int argc, char **argv) {
  // usage: test_cvole_fp <party> <port> [field]   field in {fp61, fp107, f2k}
  parse_party_and_port(argv, &party, &port);
  std::string field = (argc > 3) ? std::string(argv[3]) : "fp61";
  if (field != "fp61" && field != "fp107" && field != "f2k" && field != "z2k")
    error("unknown field: use fp61 | fp107 | f2k | z2k");

  NetIO *ios[threads];
  for (std::size_t i = 0; i < threads; ++i)
    ios[i] = new NetIO(party == ALICE ? nullptr : "127.0.0.1", port + i);

  std::cout << std::endl
            << "------------ Committed VOLE over field: " << field
            << " ------------" << std::endl
            << std::endl;

  if (field == "fp61")
    test_cvole<FP61, FP61, FP61x2>(ios, party, "fp61");
  else if (field == "fp107")
    test_cvole<FP107, FP107, FP107x2>(ios, party, "fp107");
  else if (field == "f2k")
    test_cvole<FP2x128, FP2x128, FP2x128x2>(ios, party, "f2k");
  else if (field == "z2k")
    test_cvole<Z2k64, Z2k64, Z2k64x2>(ios, party, "z2k");

  for (std::size_t i = 0; i < threads; ++i)
    delete ios[i];

  memory_usage();
  return 0;
}
