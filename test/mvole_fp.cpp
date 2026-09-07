#include "emp-tool/emp-tool.h"
#include "vole/mvole.h"
#include "vole/fields/fp61.h"
#include "vole/fields/fp61x2.h"
#include "vole/fields/fp107.h"
#include "vole/fields/fp107x2.h"
#include "vole/fields/fp2x128.h"
#include "vole/fields/z2k.h"
#include "test.h"
#include <string>
#include <vector>

using namespace emp;

// Primal-LPN n-party VOLE (vole/mvole.h) over a runtime-chosen field.
// Parties 0..n_party-1 are verifiers, party n_party is the king. Two extend
// rounds, then the Protocol-1 consistency check. Correctness is verified two
// ways: every verifier's values must equal the king's local regeneration
// (hash exchange), and the check itself pins every MAC/key relation.
//
// usage: test_mvole_fp <party> <port> [field] [cheat]
//   cheat: verifier 1 seeds its VOLE towards verifier 2 with a different seed;
//   every verifier must then abort in the check (exit code 1).
int party, port;
const std::size_t threads = 1;
const int n_party = 3;
const int rounds = 2;

template <typename FP, typename FPS>
void run(std::vector<NetIO **> &ios, int id, const char *tag, bool cheat) {
  MVoleFp<NetIO, FP, FPS> mv(id, n_party, threads, ios);
  if (cheat && id == 1) {
    mv.cheat_peer = 2;
    PRG p; p.random_block(&mv.cheat_seed, 1);
  }
  auto t0 = clock_start();
  mv.setup();
  double t_setup = time_from(t0) / 1000.0;
  std::size_t m = mv.usable;

  if (mv.is_king()) {
    std::vector<FP *> u(n_party);
    for (int i = 0; i < n_party; ++i) u[i] = new FP[m];
    for (int r = 0; r < rounds; ++r) {
      mv.extend_local(u);
      // compare with what each verifier got (its own u^i, as seen by peer 0/1)
      for (int i = 0; i < n_party; ++i) {
        block mine = Hash::hash_for_block(u[i], m * sizeof(FP)), theirs;
        ios[i][0]->recv_data(&theirs, sizeof(block));
        if (!cmpBlock(&mine, &theirs, 1)) {
          fprintf(stderr, "[%s] king: party %d round %d values differ from local regeneration\n", tag, i, r);
          std::exit(1);
        }
      }
    }
    printf("[%s] king: local regeneration matches every verifier (%d rounds, %zu values each)\n", tag, rounds, m);
    for (int i = 0; i < n_party; ++i) delete[] u[i];
    return;
  }

  FP *u = new FP[m];
  std::vector<FP *> mac(n_party, nullptr), key(n_party, nullptr);
  for (int j = 0; j < n_party; ++j)
    if (j != id) { mac[j] = new FP[m]; key[j] = new FP[m]; }
  for (int r = 0; r < rounds; ++r) {
    mv.extend(u, mac, key);
    block h = Hash::hash_for_block(u, m * sizeof(FP));
    ios[n_party][0]->send_data(&h, sizeof(block));
    ios[n_party][0]->flush();
    if (!cheat) {
      // all forward VOLEs of this party must carry the same values
      int j0 = (id == 0) ? 1 : 0;
      for (int j = 0; j < n_party; ++j) {
        if (j == id || j == j0) continue;
        for (std::size_t l = 0; l < m; ++l)
          if (!(mv.fwd_buf[j][l].getHigh() == mv.fwd_buf[j0][l].getHigh())) {
            fprintf(stderr, "[%s] party %d: values towards peers %d and %d differ\n", tag, id, j0, j);
            std::exit(1);
          }
      }
    }
  }
  mv.check();
  printf("[%s] party %d: setup %.1f ms, mesh %.1f ms, fold %.1f ms, check %.1f ms; consistency check passed\n",
         tag, id, t_setup, mv.t_mesh / 1000.0, mv.t_fold / 1000.0, mv.t_check / 1000.0);
  delete[] u;
  for (int j = 0; j < n_party; ++j) { delete[] mac[j]; delete[] key[j]; }
}

int main(int argc, char **argv) {
  parse_party_and_port(argv, &party, &port);
  std::string field = (argc > 3) ? std::string(argv[3]) : "fp61";
  bool cheat = (argc > 4) && std::string(argv[4]) == "cheat";
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
  if (party == 0)
    std::cout << std::endl << "------------ primal n-party VOLE (n=" << n_party
              << ") over field: " << field << (cheat ? " [cheat]" : "") << " ------------"
              << std::endl << std::endl;

  if (field == "fp61")       run<FP61, FP61x2>(ios, party, "fp61", cheat);
  else if (field == "fp107") run<FP107, FP107x2>(ios, party, "fp107", cheat);
  else if (field == "f2k")   run<FP2x128, FP2x128x2>(ios, party, "f2k", cheat);
  else if (field == "z2k")   run<Z2k64, Z2k64x2>(ios, party, "z2k", cheat);

  for (int j = 0; j < P; ++j) {
    if (j == party) continue;
    for (std::size_t i = 0; i < num_io; ++i) delete ios[j][i];
    delete[] ios[j];
  }
  memory_usage();
  return 0;
}
