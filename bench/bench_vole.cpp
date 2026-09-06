// Runtime benchmark for the primal silent VOLE (VoleTriple).
//
// Measures setup + in-place extend throughput at a configurable thread count.
// (Primal VOLE has a single expansion pipeline, so there is no A/B split.)
//
// usage: bench_vole <party> <port> [field] [threads] [batches]
//        field   in {fp61, fp107, f2k, z2k}   (default fp61)
//        threads                               (default 1)
//        batches = #(ot_limit-sized) extends   (default 5)
#include "emp-tool/emp-tool.h"
#include "vole/vole-fp.h"
#include "vole/fields/field_config.h"
#include <string>

using namespace emp;
using namespace std;

int party, port;

template <typename FQ, typename FP, typename FPS>
void bench_vole(NetIO **ios, std::size_t threads, std::size_t batches,
                const char *tag) {
  VoleTriple<NetIO, FP, FPS> vt(party, threads, ios);

  PRG prg;
  if (party == ALICE) {
    FP Delta;
    Delta.rand(prg);
    if (FP::PR_num_pack > 1) Delta.setHigh(FP(0));
    restrict_delta(Delta);  // ring: clamp Delta (no-op for fields)
    vt.template setup<FQ>(Delta);
  } else {
    block prog_seed;
    prg.random_block(&prog_seed, 1);
    vt.setup_prog(prog_seed);
    vt.template setup<FQ>();
  }

  std::size_t need = batches * vt.ot_limit;  // total correlations
  std::size_t mem = vt.byte_memory_need_inplace(need);

  auto start = clock_start();
  if (party == ALICE) {
    FP *buf = new FP[mem];
    vt.extend_inplace_send(buf, mem);
    delete[] buf;
  } else {
    FPS *buf = new FPS[mem];
    vt.extend_inplace_recv(buf, mem);
    delete[] buf;
  }
  double total = time_from(start) / 1000.0;  // ms

  if (party != ALICE) return;
  printf("[%s] threads=%zu batches=%zu corr=%zu\n", tag, threads, batches, need);
  printf("  extend         %9.3f ms   (%.4f us/corr, %.3f ms/batch)\n", total,
         total * 1000.0 / (double)need, total / (double)batches);
}

int main(int argc, char **argv) {
  parse_party_and_port(argv, &party, &port);
  std::string field = (argc > 3) ? std::string(argv[3]) : "fp61";
  std::size_t threads = (argc > 4) ? std::strtoull(argv[4], nullptr, 0) : 1;
  std::size_t batches = (argc > 5) ? std::strtoull(argv[5], nullptr, 0) : 5;
  if (field != "fp61" && field != "fp107" && field != "f2k" && field != "z2k")
    error("unknown field: use fp61 | fp107 | f2k | z2k");
  if (threads < 1) threads = 1;

  NetIO **ios = new NetIO *[threads];
  for (std::size_t i = 0; i < threads; ++i)
    ios[i] = new NetIO(party == ALICE ? nullptr : "127.0.0.1", port + (int)i);

  if (field == "fp61")
    bench_vole<FP61, FP61, FP61x2>(ios, threads, batches, "fp61");
  else if (field == "fp107")
    bench_vole<FP107, FP107, FP107x2>(ios, threads, batches, "fp107");
  else if (field == "f2k")
    bench_vole<FP2x128, FP2x128, FP2x128x2>(ios, threads, batches, "f2k");
  else if (field == "z2k")
    bench_vole<Z2k64, Z2k64, Z2k64x2>(ios, threads, batches, "z2k");

  for (std::size_t i = 0; i < threads; ++i) delete ios[i];
  delete[] ios;
  return 0;
}
