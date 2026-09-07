// Benchmark of the primal-LPN n-party VOLE (vole/mvole.h): n verifiers plus
// the king, all on localhost. Each verifier reports the pairwise mesh time,
// the fold time, and the consistency-check time per extend round, and the
// cost per n-party correlation (one u^i entry authenticated towards every
// other verifier).
//
// usage: bench_mvole <party> <port> [field] [threads] [n_party] [rounds]
//   party in 0..n_party (n_party = king); the mesh needs (n_party+1)^2 * io
//   ports above <port>.
#include "emp-tool/emp-tool.h"
#include "vole/mvole.h"
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
           int rounds, const char *tag) {
  MVoleFp<NetIO, FP, FPS> mv(id, n_party, threads, ios);
  auto t0 = clock_start();
  mv.setup();
  double t_setup = time_from(t0) / 1000.0;
  std::size_t m = mv.usable;

  if (mv.is_king()) {
    std::vector<FP *> u(n_party);
    for (int i = 0; i < n_party; ++i) u[i] = new FP[m];
    t0 = clock_start();
    for (int r = 0; r < rounds; ++r) mv.extend_local(u);
    double t_ext = time_from(t0) / 1000.0;
    printf("[%s] king: setup %.1f ms, local regeneration %.1f ms/round for %d verifiers\n",
           tag, t_setup, t_ext / rounds, n_party);
    for (int i = 0; i < n_party; ++i) delete[] u[i];
    return;
  }

  FP *u = new FP[m];
  std::vector<FP *> mac(n_party, nullptr), key(n_party, nullptr);
  for (int j = 0; j < n_party; ++j)
    if (j != id) { mac[j] = new FP[m]; key[j] = new FP[m]; }
  t0 = clock_start();
  for (int r = 0; r < rounds; ++r) mv.extend(u, mac, key);
  mv.check();
  double total = time_from(t0) / 1000.0;
  printf("[%s] party %d (n=%d, threads=%zu): setup %.1f ms | per round: mesh %.1f ms, fold %.1f ms, "
         "check %.1f ms (once per %d rounds) | %.4f us per n-party correlation, %zu per round\n",
         tag, id, n_party, threads, t_setup, mv.t_mesh / 1000.0 / rounds, mv.t_fold / 1000.0 / rounds,
         mv.t_check / 1000.0, rounds, total * 1000.0 / ((double)m * rounds), m);
  delete[] u;
  for (int j = 0; j < n_party; ++j) { delete[] mac[j]; delete[] key[j]; }
}

int main(int argc, char **argv) {
  int party, port;
  parse_party_and_port(argv, &party, &port);
  std::string field = (argc > 3) ? std::string(argv[3]) : "fp61";
  std::size_t threads = (argc > 4) ? (std::size_t)atoi(argv[4]) : 1;
  int n_party = (argc > 5) ? atoi(argv[5]) : 3;
  int rounds = (argc > 6) ? atoi(argv[6]) : 3;
  if (field != "fp61" && field != "fp107" && field != "f2k" && field != "z2k")
    error("unknown field: use fp61 | fp107 | f2k | z2k");

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

  if (field == "fp61")       bench<FP61, FP61x2>(ios, party, n_party, threads, rounds, "fp61");
  else if (field == "fp107") bench<FP107, FP107x2>(ios, party, n_party, threads, rounds, "fp107");
  else if (field == "f2k")   bench<FP2x128, FP2x128x2>(ios, party, n_party, threads, rounds, "f2k");
  else if (field == "z2k")   bench<Z2k64, Z2k64x2>(ios, party, n_party, threads, rounds, "z2k");

  for (int j = 0; j < P; ++j) {
    if (j == party) continue;
    for (std::size_t i = 0; i < num_io; ++i) delete ios[j][i];
    delete[] ios[j];
  }
  return 0;
}
