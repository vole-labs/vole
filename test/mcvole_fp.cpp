#include "emp-tool/emp-tool.h"
#include "vole/mcvole.h"
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

// n-party pairwise committed VOLE (paper Sec. 6.2) over a runtime-chosen field.
// Parties 0..n_party-1 are the VOLE parties; party n_party is P_king. The king
// seeds every party, publishes their commitments, and every pairwise verifier
// checks its key share against the published commitment (consistency_check_pub).
// Success = every party completes the mesh without a check failure.
const int n_party = 3;
const std::size_t target_size = 1000;  // <= ot_limit -> single round

template <typename FQ, typename FP, typename FPS>
void test_mcvole(std::vector<NetIO **> &ios, int id, const char *tag) {
  MCVoleFp<NetIO, FP, FPS> m(id, n_party, threads, ios, target_size);
  auto start = clock_start();
  m.run();
  std::cout << "[" << tag << "] party " << id << (id == n_party ? " (king)" : "")
            << " done in " << time_from(start) / 1000.0 << " ms" << std::endl;
  if (id != n_party)
    std::cout << "[" << tag << "] party " << id
              << " published-commitment checks passed" << std::endl;
}

int main(int argc, char **argv) {
  // usage: test_mcvole_fp <party> <port> [field]   party in 0..n_party (king=n_party)
  parse_party_and_port(argv, &party, &port);
  std::string field = (argc > 3) ? std::string(argv[3]) : "fp61";
  if (field != "fp61" && field != "fp107" && field != "f2k" && field != "z2k")
    error("unknown field: use fp61 | fp107 | f2k | z2k");

  const int P = n_party + 1;
  std::size_t num_io = (threads == 1) ? 2 : threads;

  // Full mesh of channels: each unordered pair gets a unique block of ports.
  // Lower-id party is the server, higher-id party is the client.
  std::vector<NetIO **> ios(P, nullptr);
  for (int j = 0; j < P; ++j) {
    if (j == party) continue;
    ios[j] = new NetIO *[num_io];
    int lo = std::min(party, j), hi = std::max(party, j);
    for (std::size_t i = 0; i < num_io; ++i) {
      int port_ = port + (lo * P + hi) * (int)num_io + (int)i;
      if (party < j)  // this party = server (lower id)
        ios[j][i] = new NetIO(nullptr, port_);
      else            // this party = client (higher id)
        ios[j][i] = new NetIO("127.0.0.1", port_);
    }
  }

  if (party == 0)
    std::cout << std::endl
              << "------------ n-party committed VOLE (n=" << n_party
              << ") over field: " << field << " ------------" << std::endl
              << std::endl;

  if (field == "fp61")
    test_mcvole<FP61, FP61, FP61x2>(ios, party, "fp61");
  else if (field == "fp107")
    test_mcvole<FP107, FP107, FP107x2>(ios, party, "fp107");
  else if (field == "f2k")
    test_mcvole<FP2x128, FP2x128, FP2x128x2>(ios, party, "f2k");
  else if (field == "z2k")
    test_mcvole<Z2k64, Z2k64, Z2k64x2>(ios, party, "z2k");

  for (int j = 0; j < P; ++j) {
    if (j == party) continue;
    for (std::size_t i = 0; i < num_io; ++i) delete ios[j][i];
    delete[] ios[j];
  }

  memory_usage();
  return 0;
}
