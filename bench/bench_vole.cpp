// Runtime benchmark for the primal silent VOLE (VoleTriple).
//
// Measures setup + in-place extend throughput at a configurable thread count.
// (Primal VOLE has a single expansion pipeline, so there is no A/B split.)
//
// usage: bench_vole <party> <port> [field] [threads] [batches] [params]
//        field   in {fp61, fp107, f2k, z2k, f2}  (default fp61; f2 = Ferret COT,
//                  values in F_2, always uses fp_ferret_f2)
//        threads                               (default 1)
//        batches = #(ot_limit-sized) extends   (default 5)
//        params  in {wolverine, ferret_b13}    (default wolverine)
#include "emp-tool/emp-tool.h"
#include "vole/vole-fp.h"
#include "vole/fields/field_config.h"
#include <string>

using namespace emp;
using namespace std;

int party, port;

template <typename FQ, typename FP, typename FPS>
void bench_vole(NetIO **ios, std::size_t threads, std::size_t batches,
                const char *tag, const PrimalLPNParameterFp61 &param) {
  VoleTriple<NetIO, FP, FPS> vt(party, threads, ios, param);

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

  // Allocate (and touch) the output buffer outside the timed region: for
  // batches * ot_limit outputs this is hundreds of MB of single-threaded
  // construction and page faults, which would otherwise dominate the
  // multi-threaded numbers (emp-ot's benches allocate before timing too).
  FP *buf_k = nullptr; FPS *buf_vm = nullptr;
  if (party == ALICE) buf_k = new FP[mem]; else buf_vm = new FPS[mem];
  auto start = clock_start();
  if (party == ALICE) vt.extend_inplace_send(buf_k, mem);
  else vt.extend_inplace_recv(buf_vm, mem);
  double total = time_from(start) / 1000.0;  // ms
  delete[] buf_k; delete[] buf_vm;

  printf("[%s %s] threads=%zu batches=%zu corr=%zu (n=%zu t=%zu k=%zu)\n", tag,
         party == ALICE ? "ALICE" : "BOB", threads, batches, need, param.n, param.t, param.k);
  printf("  extend         %9.3f ms   (%.4f us/corr, %.3f ms/batch)\n", total,
         total * 1000.0 / (double)need, total / (double)batches);
  printf("  per batch: cot %.1f  mpfss %.1f  lpn %.1f  copy %.1f  | extend call %.1f ms\n",
         vt.t_cot / 1000.0 / batches, vt.t_mpfss / 1000.0 / batches,
         vt.t_lpn / 1000.0 / batches, vt.t_copy / 1000.0 / batches, vt.t_round / 1000.0 / batches);
}

int main(int argc, char **argv) {
  parse_party_and_port(argv, &party, &port);
  std::string field = (argc > 3) ? std::string(argv[3]) : "fp61";
  std::size_t threads = (argc > 4) ? std::strtoull(argv[4], nullptr, 0) : 1;
  std::size_t batches = (argc > 5) ? std::strtoull(argv[5], nullptr, 0) : 5;
  std::string pname = (argc > 6) ? std::string(argv[6]) : "wolverine";
  if (pname != "wolverine" && pname != "ferret_b13")
    error("unknown params: use wolverine | ferret_b13");
  const PrimalLPNParameterFp61 &param = (pname == "ferret_b13") ? fp_ferret_b13 : fp_default;
  if (field != "fp61" && field != "fp107" && field != "f2k" && field != "z2k" && field != "f2")
    error("unknown field: use fp61 | fp107 | f2k | z2k | f2");
  if (threads < 1) threads = 1;

  NetIO **ios = new NetIO *[threads];
  for (std::size_t i = 0; i < threads; ++i)
    ios[i] = new NetIO(party == ALICE ? nullptr : "127.0.0.1", port + (int)i);

  if (field == "fp61")
    bench_vole<FP61, FP61, FP61x2>(ios, threads, batches, "fp61", param);
  else if (field == "fp107")
    bench_vole<FP107, FP107, FP107x2>(ios, threads, batches, "fp107", param);
  else if (field == "f2k")
    bench_vole<FP2x128, FP2x128, FP2x128x2>(ios, threads, batches, "f2k", param);
  else if (field == "z2k")
    bench_vole<Z2k64, Z2k64, Z2k64x2>(ios, threads, batches, "z2k", param);
  else if (field == "f2")
    bench_vole<F2kKey, F2kKey, F2Auth>(ios, threads, batches, "f2", fp_ferret_f2);

  for (std::size_t i = 0; i < threads; ++i) delete ios[i];
  delete[] ios;
  return 0;
}
