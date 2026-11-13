// Two-party check of COPE (vole/cope.h), the base correlation under the base
// sVOLE: the sender holds Delta, the receiver inputs u, and the outputs satisfy
//   recv_out[j] = send_out[j] + Delta * u[j].
// Runs fp107 (built-in integer storage, the original COPE code) and p192 (limb
// storage, CopeFp::limb_extend_*), including edge inputs u in {0, 1, p-1} and
// a Delta with both the low and the top bit set (every bit of the
// decomposition is exercised).
//
// usage: test_cope <party> <port> [field]   party 1 = sender (ALICE), 2 = receiver
#include "emp-tool/emp-tool.h"
#include "vole/cope.h"
#include "vole/fields/field_config.h"
#include "test.h"
#include <string>

using namespace emp;

template <typename FP>
bool bit_of(const FP &x, std::size_t i) {
  if constexpr (std::is_same_v<decltype(FP::val), U192>) return x.val.bit(i);
  else return (x.val >> i) & 1;
}

template <typename FP>
int run(NetIO *io, int party, const char *tag) {
  const std::size_t m = FP::slot_stride_bits, size = 1000;
  CopeFp<NetIO, FP> cope(party, io);
  std::vector<FP> out(size);
  PRG prg;
  if (party == ALICE) {
    // resample until bit 0 and the top bit (m - 1) of Delta are both set
    FP Delta;
    do Delta.rand(prg); while (!(bit_of(Delta, 0) && bit_of(Delta, m - 1)));
    cope.initialize(Delta, m);
    cope.extend64(out.data(), size, m, 1);
    std::vector<FP> u(size), rout(size);
    io->recv_data(u.data(), size * sizeof(FP));
    io->recv_data(rout.data(), size * sizeof(FP));
    std::size_t bad = 0;
    for (std::size_t j = 0; j < size; ++j)
      if (rout[j] != out[j] + Delta * u[j]) ++bad;
    printf("[%s] COPE relation: %zu / %zu mismatches\n", tag, bad, size);
    return bad == 0 ? 0 : 1;
  } else {
    std::vector<FP> u(size);
    for (std::size_t j = 0; j < size; ++j) u[j].rand(prg);
    u[0] = FP(0);
    u[1] = FP(1);
    u[2] = FP(-1);  // p - 1
    cope.initialize(m);
    cope.extend64(out.data(), u.data(), size, m, 1);
    io->send_data(u.data(), size * sizeof(FP));
    io->send_data(out.data(), size * sizeof(FP));
    io->flush();
    return 0;
  }
}

int main(int argc, char **argv) {
  int party, port;
  parse_party_and_port(argv, &party, &port);
  std::string field = (argc > 3) ? std::string(argv[3]) : "p192";
  NetIO *io = new NetIO(party == ALICE ? nullptr : "127.0.0.1", port);
  int r;
  if (field == "fp107") r = run<FP107>(io, party, "fp107");
  else if (field == "p192") r = run<FP192>(io, party, "p192");
  else error("unknown field: use fp107 | p192");
  delete io;
  return r;
}
