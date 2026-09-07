// Benchmark of the dual-LPN n-party committed VOLE (vole/mcvole.h), the
// counterpart of bench_mvole: n verifiers plus the king on localhost. Each
// verifier reports its pairwise mesh time (all committed VOLEs with every
// peer, including the published-commitment checks) and the cost per n-party
// correlation, excluding setup, so the numbers compare directly with
// bench_mvole's mesh + fold + check.
//
// usage: bench_mcvole <party> <port> [field] [threads] [n_party] [log2_target]
#include "emp-tool/emp-tool.h"
#include "vole/mcvole.h"
#include "vole/fields/fp61.h"
#include "vole/fields/fp61x2.h"
#include "vole/fields/fp107.h"
#include "vole/fields/fp107x2.h"
#include "vole/fields/fp2x128.h"
#include "vole/fields/z2k.h"
#include "test/test.h"
#include <string>
#include <vector>

using namespace emp;

template <typename FP, typename FPS>
void bench(std::vector<NetIO **> &ios, int id, int n_party, std::size_t threads,
           std::size_t target, const char *tag) {
  MCVoleFp<NetIO, FP, FPS> mc(id, n_party, threads, ios, target);
  std::size_t M = (mc.param.t + 1) + (mc.param.t_com + 1);
  std::size_t per_round = mc.param.n - M;
  std::size_t outputs = mc.rounds * per_round;
  auto t0 = clock_start();
  if (mc.is_king()) {
    mc.run_king();
    printf("[%s] king: seeds + local regeneration + publish %.1f ms for %d verifiers (%zu rounds)\n",
           tag, time_from(t0) / 1000.0, n_party, mc.rounds);
    return;
  }
  // run_party() with a timer around the extend mesh.
  ios[n_party][0]->recv_data(&mc.prog_seed, sizeof(block));
  mc.published_com.resize(n_party);
  for (int i = 0; i < n_party; ++i) {
    mc.published_com[i].resize(mc.cn);
    ios[n_party][0]->recv_data(mc.published_com[i].data(), mc.cn * sizeof(FP));
  }
  mc.voleForward.resize(n_party, nullptr);
  mc.voleInverse.resize(n_party, nullptr);
  mc.mesh_run(/*is_setup=*/true);
  double t_setup = time_from(t0) / 1000.0;
  t0 = clock_start();
  mc.mesh_run(/*is_setup=*/false);
  double t_mesh = time_from(t0) / 1000.0;
  printf("[%s] party %d (n=%d, threads=%zu): setup %.1f ms | mesh %.1f ms for %zu rounds "
         "(%.1f ms/round) | %.4f us per n-party correlation, %zu per round\n",
         tag, id, n_party, threads, t_setup, t_mesh, mc.rounds, t_mesh / mc.rounds,
         t_mesh * 1000.0 / (double)outputs, per_round);
}

int main(int argc, char **argv) {
  int party, port;
  parse_party_and_port(argv, &party, &port);
  std::string field = (argc > 3) ? std::string(argv[3]) : "fp61";
  std::size_t threads = (argc > 4) ? (std::size_t)atoi(argv[4]) : 1;
  int n_party = (argc > 5) ? atoi(argv[5]) : 3;
  int log2_target = (argc > 6) ? atoi(argv[6]) : 20;
  if (field != "fp61" && field != "fp107" && field != "f2k" && field != "z2k")
    error("unknown field: use fp61 | fp107 | f2k | z2k");
  std::size_t target = 1ull << log2_target;

  const int P = n_party + 1;
  std::size_t num_io = (threads == 1) ? 2 : threads;
  std::vector<NetIO **> ios(P, nullptr);
  for (int j = 0; j < P; ++j) {
    if (j == party) continue;
    ios[j] = new NetIO *[num_io];
    int lo = std::min(party, j), hi = std::max(party, j);
    for (std::size_t i = 0; i < num_io; ++i) {
      int port_ = port + (lo * P + hi) * (int)num_io + (int)i;
      if (party < j) ios[j][i] = new NetIO(nullptr, port_);
      else           ios[j][i] = new NetIO("127.0.0.1", port_);
    }
  }

  if (field == "fp61")       bench<FP61, FP61x2>(ios, party, n_party, threads, target, "fp61");
  else if (field == "fp107") bench<FP107, FP107x2>(ios, party, n_party, threads, target, "fp107");
  else if (field == "f2k")   bench<FP2x128, FP2x128x2>(ios, party, n_party, threads, target, "f2k");
  else if (field == "z2k")   bench<Z2k64, Z2k64x2>(ios, party, n_party, threads, target, "z2k");

  for (int j = 0; j < P; ++j) {
    if (j == party) continue;
    for (std::size_t i = 0; i < num_io; ++i) delete ios[j][i];
    delete[] ios[j];
  }
  return 0;
}
