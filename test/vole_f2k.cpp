// Two-party test of the Ferret-style VOLE interface (vole/vole_f2k.h):
//   * F2kVole<NetIO> over GF(2^128) with block arguments: random VOLE served
//     in odd-sized requests that straddle round refills, then chosen-input
//     VOLE; both checked against mac == key ^ gfmul(val, Delta);
//   * the generic RVole over fp61, same checks with field arithmetic.
// usage: test_vole_f2k <party> <port>
#include "emp-tool/emp-tool.h"
#include "vole/vole_f2k.h"
#include "vole/fields/fp61.h"
#include "vole/fields/fp61x2.h"
#include <vector>
#include "test.h"

using namespace emp;
int party, port;
const std::size_t threads = 1;

// Sizes that cross a refill boundary: ot_limit - 5, then 1000, then 7.
static std::vector<std::size_t> sizes(std::size_t ot_limit) {
  return {ot_limit - 5, 1000, 7, ot_limit + 3};
}

// ALICE ships Delta and keys; BOB verifies mac == key ^ gfmul(val, Delta).
static void check_f2k(NetIO *io, const block *key, const block *val, const block *mac,
                      block Delta, std::size_t n, const char *what) {
  if (party == ALICE) {
    io->send_data(&Delta, sizeof(block));
    io->send_data(key, n * sizeof(block));
    io->flush();
  } else {
    block d; std::vector<block> k(n);
    io->recv_data(&d, sizeof(block));
    io->recv_data(k.data(), n * sizeof(block));
    for (std::size_t i = 0; i < n; ++i) {
      block t; gfmul(val[i], d, &t); t = t ^ k[i];
      if (!cmpBlock(&t, &mac[i], 1)) error(what);
    }
    printf("  [f2k] %s: %zu correlations ok\n", what, n);
  }
}

void test_f2k(NetIO **ios) {
  F2kVole<NetIO> v(party, threads, ios);
  v.setup();
  std::vector<std::size_t> sz = sizes(v.ot_limit());
  std::size_t total = 0; for (auto s : sz) total += s;
  std::vector<block> key(total), val(total), mac(total);
  std::size_t off = 0;
  for (auto s : sz) {
    if (party == ALICE) v.rvole(key.data() + off, s);
    else v.rvole(val.data() + off, mac.data() + off, s);
    off += s;
  }
  check_f2k(ios[0], key.data(), val.data(), mac.data(), v.delta(), total, "random VOLE across refills");

  // chosen input
  std::size_t n = 12345;
  std::vector<block> x(n), k2(n), m2(n);
  PRG prg; prg.random_block(x.data(), n);
  if (party == ALICE) v.vole_send(k2.data(), n);
  else v.vole_recv(x.data(), m2.data(), n);
  check_f2k(ios[0], k2.data(), x.data(), m2.data(), v.delta(), n, "chosen-input VOLE");
}

void test_fp61(NetIO **ios) {
  RVole<NetIO, FP61, FP61x2> v(party, threads, ios);
  PRG prg;
  if (party == ALICE) { FP61 d; d.rand(prg); v.setup(d); } else v.setup();
  std::vector<std::size_t> sz = sizes(v.ot_limit());
  std::size_t total = 0; for (auto s : sz) total += s;
  std::vector<FP61> key(total); std::vector<FP61x2> vm(total);
  std::size_t off = 0;
  for (auto s : sz) {
    if (party == ALICE) v.rvole_send(key.data() + off, s); else v.rvole_recv(vm.data() + off, s);
    off += s;
  }
  std::size_t n = 12345;
  std::vector<FP61> x(n), k2(n), m2(n);
  for (auto &e : x) e.rand(prg);
  if (party == ALICE) v.vole_send(k2.data(), n); else v.vole_recv(x.data(), m2.data(), n);

  if (party == ALICE) {
    ios[0]->send_data(&v.Delta, sizeof(FP61));
    ios[0]->send_data(key.data(), total * sizeof(FP61));
    ios[0]->send_data(k2.data(), n * sizeof(FP61));
    ios[0]->flush();
  } else {
    FP61 d; std::vector<FP61> k(total), kk(n);
    ios[0]->recv_data(&d, sizeof(FP61));
    ios[0]->recv_data(k.data(), total * sizeof(FP61));
    ios[0]->recv_data(kk.data(), n * sizeof(FP61));
    for (std::size_t i = 0; i < total; ++i)
      if (!(vm[i].getLow() == k[i] + vm[i].getHigh() * d)) error("fp61 random VOLE");
    for (std::size_t i = 0; i < n; ++i)
      if (!(m2[i] == kk[i] + x[i] * d)) error("fp61 chosen-input VOLE");
    printf("  [fp61] random (%zu) and chosen-input (%zu) VOLE ok\n", total, n);
  }
}

int main(int argc, char **argv) {
  parse_party_and_port(argv, &party, &port);
  NetIO *ios[threads];
  for (std::size_t i = 0; i < threads; ++i)
    ios[i] = new NetIO(party == ALICE ? nullptr : "127.0.0.1", port + i);
  test_f2k(ios);
  test_fp61(ios);
  for (std::size_t i = 0; i < threads; ++i) delete ios[i];
  if (party == BOB) printf("Ferret-style VOLE interface: all checks passed\n");
  return 0;
}
