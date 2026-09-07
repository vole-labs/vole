// Runtime benchmark for the committed VOLE (CVoleFp).
//
// Unlike test/cvole_fp.cpp (which checks correctness with a tiny param), this
// measures throughput at a configurable scale and thread count, and breaks the
// extend time into the two pipelines discussed in the design:
//   A (VOLE)       : COT + MPFSS(e_u) + LPN H      -> x
//   B (commitment) : COT + MPFSS(e_r) + [H_u|H_r]*[e_u||e_r] -> com
//   shared         : the accumulate over e_u (feeds H)
// (The committed-VOLE opening proof, consistency_check_*, is a separate step and
// is NOT part of extend, so it is not measured here.)
//
// usage: bench_cvole <party> <port> [field] [threads] [log2_target]
//        field   in {fp61, fp107, f2k, z2k}   (default fp61)
//        threads                               (default 1)
//        log2_target = log2(#correlations)     (default 20)
#include "emp-tool/emp-tool.h"
#include "vole/cvole.h"
#include "vole/fields/field_config.h"
#include <string>

using namespace emp;
using namespace std;

int party, port;

template <typename FQ, typename FP, typename FPS>
void bench_cvole(NetIO **ios, std::size_t threads, std::size_t target,
                 const char *tag) {
  CVoleFpParam param = cvole_fp_param_for(target);
  CVoleFp<NetIO, FP, FPS> cvole(party, threads, ios, param);

  PRG prg;
  if (party == ALICE) {
    FP Delta;
    Delta.rand(prg);
    if (FP::PR_num_pack > 1) Delta.setHigh(FP(0));
    restrict_delta(Delta);  // ring: clamp Delta (no-op for fields)
    cvole.setup(Delta);
  } else {
    block prog_seed;
    prg.random_block(&prog_seed, 1);
    cvole.setup_prog(prog_seed);
    cvole.setup();
  }

  std::size_t rounds = cvole.rounds_for(target);
  std::size_t xn = cvole.x_buf_size(rounds);
  std::size_t cn = cvole.com_buf_size(rounds);
  std::size_t corr = cvole.x_usable(rounds);   // usable VOLE correlations

  // Buffers are allocated (and touched) outside the timed region, as in
  // emp-ot's benches; their single-threaded construction is not protocol work.
  FP *xk = nullptr, *ck = nullptr; FPS *xv = nullptr, *cv = nullptr;
  if (party == ALICE) { xk = new FP[xn]; ck = new FP[cn]; }
  else { xv = new FPS[xn]; cv = new FPS[cn]; }
  auto start = clock_start();
  if (party == ALICE) cvole.extend_inplace_send(xk, ck, rounds);
  else cvole.extend_inplace_recv(xv, cv, rounds);
  double total = time_from(start) / 1000.0;  // ms
  delete[] xk; delete[] ck; delete[] xv; delete[] cv;

  // Only one party needs to print; ALICE here. Both see equal phase splits since
  // the parties are synchronized inside each MPFSS.
  if (party != ALICE) return;
  printf("[%s] threads=%zu rounds=%zu corr=%zu\n", tag, threads, rounds, corr);
  printf("  total extend   %9.3f ms   (%.4f us/corr)\n", total,
         total * 1000.0 / (double)corr);
  printf("  A  VOLE        %9.3f ms\n", cvole.t_vole / 1000.0);
  printf("  B  commitment  %9.3f ms\n", cvole.t_com / 1000.0);
  printf("  shared accum   %9.3f ms\n", cvole.t_accum / 1000.0);
}

int main(int argc, char **argv) {
  parse_party_and_port(argv, &party, &port);
  std::string field = (argc > 3) ? std::string(argv[3]) : "fp61";
  std::size_t threads = (argc > 4) ? std::strtoull(argv[4], nullptr, 0) : 1;
  std::size_t log2_target = (argc > 5) ? std::strtoull(argv[5], nullptr, 0) : 20;
  std::size_t target = 1ull << log2_target;
  if (field != "fp61" && field != "fp107" && field != "f2k" && field != "z2k")
    error("unknown field: use fp61 | fp107 | f2k | z2k");
  if (threads < 1) threads = 1;

  NetIO **ios = new NetIO *[threads];
  for (std::size_t i = 0; i < threads; ++i)
    ios[i] = new NetIO(party == ALICE ? nullptr : "127.0.0.1", port + (int)i);

  if (field == "fp61")
    bench_cvole<FP61, FP61, FP61x2>(ios, threads, target, "fp61");
  else if (field == "fp107")
    bench_cvole<FP107, FP107, FP107x2>(ios, threads, target, "fp107");
  else if (field == "f2k")
    bench_cvole<FP2x128, FP2x128, FP2x128x2>(ios, threads, target, "f2k");
  else if (field == "z2k")
    bench_cvole<Z2k64, Z2k64, Z2k64x2>(ios, threads, target, "z2k");

  for (std::size_t i = 0; i < threads; ++i) delete ios[i];
  delete[] ios;
  return 0;
}
