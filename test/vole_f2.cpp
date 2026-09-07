// Two-party test of the F_2-value instantiation (Ferret correlated OT) built on
// VoleTriple<NetIO, F2kKey, F2Auth> through the FerretVole interface:
//   * rcot in odd sizes across refills: every receiver block R satisfies
//     R == K xor (LSB(R) ? Delta : 0), sender keys have LSB 0, Delta has LSB 1;
//   * recv_cot with chosen bits: R == K' xor (b ? Delta : 0).
// usage: test_vole_f2 <party> <port>
#include "emp-tool/emp-tool.h"
#include "vole/vole_f2k.h"
#include <vector>

using namespace emp;
int party, port;
const std::size_t threads = 1;

// `b == nullptr` means random COT: then the receiver's bit is LSB(R) and the
// sender's keys must have LSB 0. After the chosen-bit online phase the sender's
// corrected keys K' = K xor d*Delta carry LSB d (as in Ferret), so that
// invariant is only checked for the random phase.
static void check(NetIO *io, const block *key, const block *r, const bool *b, block Delta,
                  std::size_t n, const char *what) {
  if (party == ALICE) {
    if (b == nullptr)
      for (std::size_t i = 0; i < n; ++i) if (getLSB(key[i])) error("sender key with LSB 1");
    if (!getLSB(Delta)) error("Delta with LSB 0");
    io->send_data(&Delta, sizeof(block));
    io->send_data(key, n * sizeof(block));
    io->flush();
  } else {
    block d; std::vector<block> k(n);
    io->recv_data(&d, sizeof(block));
    io->recv_data(k.data(), n * sizeof(block));
    std::size_t ones = 0;
    for (std::size_t i = 0; i < n; ++i) {
      bool bit = b ? b[i] : getLSB(r[i]);
      block t = bit ? (k[i] ^ d) : k[i];
      if (!cmpBlock(&t, &r[i], 1)) error(what);
      ones += bit;
    }
    printf("  [f2] %s: %zu COTs ok (%.3f ones)\n", what, n, (double)ones / n);
  }
}

int main(int argc, char **argv) {
  parse_party_and_port(argv, &party, &port);
  NetIO *ios[threads];
  for (std::size_t i = 0; i < threads; ++i)
    ios[i] = new NetIO(party == ALICE ? nullptr : "127.0.0.1", port + i);

  FerretVole<NetIO> v(party, threads, ios);
  auto t0 = clock_start();
  v.setup();
  if (party == BOB) printf("  [f2] setup %.1f ms (ot_limit %zu)\n", time_from(t0) / 1000.0, v.ot_limit());

  std::vector<std::size_t> sz = {v.ot_limit() - 5, 1000, 7, v.ot_limit() + 3};
  std::size_t total = 0; for (auto s : sz) total += s;
  std::vector<block> buf(total);
  std::size_t off = 0;
  t0 = clock_start();
  for (auto s : sz) { v.rcot(buf.data() + off, s); off += s; }
  if (party == BOB) printf("  [f2] rcot %zu in %.1f ms\n", total, time_from(t0) / 1000.0);
  check(ios[0], buf.data(), buf.data(), nullptr, v.delta(), total, "random COT across refills");

  std::size_t n = 12345;
  std::vector<block> c(n);
  bool *b = new bool[n];
  PRG prg; prg.random_bool(b, n);
  if (party == ALICE) v.send_cot(c.data(), n); else v.recv_cot(c.data(), b, n);
  check(ios[0], c.data(), c.data(), b, v.delta(), n, "chosen-bit COT");
  delete[] b;

  for (std::size_t i = 0; i < threads; ++i) delete ios[i];
  if (party == BOB) printf("F_2 VOLE (Ferret COT) checks passed\n");
  return 0;
}
