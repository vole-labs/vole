#include "emp-tool/emp-tool.h"
#include "vole/ncvole.h"
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
const int n_client = 2;

// Target number of committed VOLE outputs; overridable via argv[4] for
// benchmarking (default = 2^20 to amortise setup; 2^14 is fast for smoke tests).
std::size_t target_size = 1ull << 20;

// Server (committer) runs committed VOLE with every client. The committed value
// is materialized at the first client and reused for the rest; each client's
// own consistency check (client side) enforces it uses that one published
// commitment. Buffers are managed inside the driver.
template <typename FQ, typename FP, typename FPS>
void run_server(std::vector<NetIO **> &ios, const char *tag) {
  ProgNCVoleFp<NetIO, FP, FPS> nvole(0, n_client, threads, target_size);
  // usable correlations per server-client extend (rounds * ot_limit)
  std::size_t M = (nvole.param.t + 1) + (nvole.param.t_com + 1);
  std::size_t ot_limit = nvole.param.n - M;
  std::size_t rounds = (target_size + ot_limit - 1) / ot_limit;
  std::size_t usable = rounds * ot_limit;
  for (int i = 0; i < n_client; ++i) {
    auto t0 = clock_start();
    nvole.extend_server(i, ios[i]);
    double ms = time_from(t0) / 1000.0;
    std::cout << "[" << tag << "] server <-> client " << i
              << " : " << ms << " ms, "
              << (ms * 1000.0 / (double)usable) << " us/corr ("
              << usable << " correlations)" << std::endl;
  }
  std::cout << "[" << tag << "] server: one committed input published to all "
            << n_client << " clients" << std::endl;
}

template <typename FQ, typename FP, typename FPS>
void run_client(NetIO **ios, const char *tag) {
  ProgNCVoleFp<NetIO, FP, FPS> nvole(party, n_client, threads, target_size);
  std::size_t M = (nvole.param.t + 1) + (nvole.param.t_com + 1);
  std::size_t ot_limit = nvole.param.n - M;
  std::size_t rounds = (target_size + ot_limit - 1) / ot_limit;
  std::size_t usable = rounds * ot_limit;
  auto t0 = clock_start();
  nvole.extend_client(ios);
  double ms = time_from(t0) / 1000.0;
  std::cout << "[" << tag << "] client " << party
            << " consistency check passed : " << ms << " ms, "
            << (ms * 1000.0 / (double)usable) << " us/corr ("
            << usable << " correlations)" << std::endl;
}

int main(int argc, char **argv) {
  // usage: server: test_ncvole_fp 0 <port> [field]
  //        client: test_ncvole_fp <1..n_client> <port> [field]
  //        field in {fp61, fp107, f2k}
  parse_party_and_port(argv, &party, &port);
  std::string field = (argc > 3) ? std::string(argv[3]) : "fp61";
  if (field != "fp61" && field != "fp107" && field != "f2k" && field != "z2k")
    error("unknown field: use fp61 | fp107 | f2k | z2k");
  if (argc > 4) target_size = std::strtoull(argv[4], nullptr, 0);

  std::cout << std::endl
            << "------------ Multi-client committed VOLE over field: " << field
            << " ------------" << std::endl
            << std::endl;

  if (party == 0) {
    std::vector<NetIO **> ios(n_client);
    for (int i = 0; i < n_client; ++i) {
      ios[i] = new NetIO *[threads];
      for (std::size_t j = 0; j < threads; ++j)
        ios[i][j] = new NetIO(nullptr, port + i * threads + j);
    }

    if (field == "fp61")       run_server<FP61, FP61, FP61x2>(ios, "fp61");
    else if (field == "fp107") run_server<FP107, FP107, FP107x2>(ios, "fp107");
    else if (field == "f2k")   run_server<FP2x128, FP2x128, FP2x128x2>(ios, "f2k");
    else if (field == "z2k")   run_server<Z2k64, Z2k64, Z2k64x2>(ios, "z2k");

    for (int i = 0; i < n_client; ++i) {
      for (std::size_t j = 0; j < threads; ++j) delete ios[i][j];
      delete[] ios[i];
    }
  } else {
    NetIO **ios = new NetIO *[threads];
    for (std::size_t j = 0; j < threads; ++j)
      ios[j] = new NetIO("127.0.0.1", port + (party - 1) * threads + j);

    if (field == "fp61")       run_client<FP61, FP61, FP61x2>(ios, "fp61");
    else if (field == "fp107") run_client<FP107, FP107, FP107x2>(ios, "fp107");
    else if (field == "f2k")   run_client<FP2x128, FP2x128, FP2x128x2>(ios, "f2k");
    else if (field == "z2k")   run_client<Z2k64, Z2k64, Z2k64x2>(ios, "z2k");

    for (std::size_t j = 0; j < threads; ++j) delete ios[j];
    delete[] ios;
  }

  memory_usage();
  return 0;
}
