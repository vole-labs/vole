#include "emp-tool/emp-tool.h"
#include "vole/vole-fp.h"
// All supported fields are compiled in; the choice is made at run time (argv[3])
// so one binary can be run over any of them.
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

// Field-generic VOLE smoke test.
// FQ : clear (single-element) field used to seed the receiver / consistency.
// FP : single-element key/MAC field (sender side).
// FPS: packed (value||MAC) bundle (receiver side).
template <typename FQ, typename FP, typename FPS>
void test_vole_triple(NetIO *ios[threads], int party, const char *tag) {
  VoleTriple<NetIO, FP, FPS> vtriple(party, threads, ios);

  FP Delta;
  block prog_seed;
  PRG prg;
  if (party == ALICE) {
    Delta.rand(prg);
    if (FP::PR_num_pack > 1)
      Delta.setHigh(FP(0));
    restrict_delta(Delta);  // ring: clamp Delta to its s-bit sub-ring (no-op for fields)

    auto start = clock_start();
    vtriple.template setup<FQ>(Delta);
    std::cout << "[" << tag << "] setup " << time_from(start) / 1000 << " ms"
              << std::endl;

    vtriple.check_triple(Delta, vtriple.pre_yz_send, vtriple.param.n_pre, ios[0]);
  } else {
    prg.random_block(&prog_seed, 1);
    auto start = clock_start();
    vtriple.setup_prog(prog_seed);
    vtriple.template setup<FQ>();
    std::cout << "[" << tag << "] setup " << time_from(start) / 1000 << " ms"
              << std::endl;

    vtriple.check_triple(vtriple.pre_yz_recv, vtriple.param.n_pre, ios[0]);
  }

  // in-place extension over several batches
  uint64_t batch_n = 5;
  uint64_t triple_need_inplace = batch_n * vtriple.ot_limit;
  uint64_t memory_need = vtriple.byte_memory_need_inplace(triple_need_inplace);
  if (party == ALICE) {
    FP *buf = new FP[memory_need];
    auto start = clock_start();
    vtriple.extend_inplace_send(buf, memory_need);
    std::cout << "[" << tag << "] extend "
              << time_from(start) / 1000 / batch_n << " ms" << std::endl;
    vtriple.check_triple(Delta, buf, memory_need, ios[0]);
    delete[] buf;
  } else {
    FPS *buf = new FPS[memory_need];
    auto start = clock_start();
    vtriple.extend_inplace_recv(buf, memory_need);
    std::cout << "[" << tag << "] extend "
              << time_from(start) / 1000 / batch_n << " ms" << std::endl;
    vtriple.check_triple(buf, memory_need, ios[0]);
    delete[] buf;
  }

  // repeated extension
  std::size_t test_n = 4;
  if (party == ALICE) {
    FP *buf = new FP[memory_need];
    for (std::size_t i = 0; i < test_n; ++i) {
      vtriple.extend_inplace_send(buf, memory_need);
      vtriple.check_triple(Delta, buf, memory_need, ios[0]);
    }
    delete[] buf;
  } else {
    FPS *buf = new FPS[memory_need];
    for (std::size_t i = 0; i < test_n; ++i) {
      vtriple.extend_inplace_recv(buf, memory_need);
      vtriple.check_triple(buf, memory_need, ios[0]);
    }
    delete[] buf;
  }
  std::cout << "[" << tag << "] pass repeated extension checks" << std::endl;
}

int main(int argc, char **argv) {
  // usage: test_vole_fp <party> <port> [field]   field in {fp61, fp107, f2k, z2k}
  parse_party_and_port(argv, &party, &port);
  std::string field = (argc > 3) ? std::string(argv[3]) : "fp61";
  if (field != "fp61" && field != "fp107" && field != "f2k" && field != "z2k")
    error("unknown field: use fp61 | fp107 | f2k | z2k");

  NetIO *ios[threads];
  for (std::size_t i = 0; i < threads; ++i)
    ios[i] = new NetIO(party == ALICE ? nullptr : "127.0.0.1", port + i);

  std::cout << std::endl
            << "------------ VOLE over field: " << field << " ------------"
            << std::endl
            << std::endl;

  if (field == "fp61")
    test_vole_triple<FP61, FP61, FP61x2>(ios, party, "fp61");
  else if (field == "fp107")
    test_vole_triple<FP107, FP107, FP107x2>(ios, party, "fp107");
  else if (field == "f2k")
    test_vole_triple<FP2x128, FP2x128, FP2x128x2>(ios, party, "f2k");
  else if (field == "z2k")
    test_vole_triple<Z2k64, Z2k64, Z2k64x2>(ios, party, "z2k");

  for (std::size_t i = 0; i < threads; ++i)
    delete ios[i];

  memory_usage();
  return 0;
}
