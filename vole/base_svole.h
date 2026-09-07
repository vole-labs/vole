#ifndef BASE_VOLE_H__
#define BASE_VOLE_H__

#include "emp-tool/emp-tool.h"
#include "vole/cope.h"
#include "vole/fields/field_config.h"

template <typename IO, typename FP> 
class BaseSvoleFp {
public:
  int party;
  std::size_t m;
  IO *io;
  CopeFp<IO, FP> *cope;
  FP Delta;
  PRG prg;
  std::size_t num_pack;

  // SENDER
  BaseSvoleFp(int party, IO *io, FP Delta) {
    this->party = party;
    this->io = io;
    cope = new CopeFp<IO, FP>(party, io);
    this->Delta = Delta;
    num_pack = (FP::PR_num_pack>1)?(FP::PR_num_pack-1):1;
    cope->initialize(Delta, FP::slot_stride_bits*num_pack);
  }

  // RECEIVER
  BaseSvoleFp(int party, IO *io) {
    this->party = party;
    this->io = io;
    cope = new CopeFp<IO, FP>(party, io);
    num_pack = (FP::PR_num_pack>1)?(FP::PR_num_pack-1):1;
    cope->initialize(FP::slot_stride_bits*num_pack);
  }

  ~BaseSvoleFp() { delete cope; }

  // sender
  void compute_send64(FP *share, std::size_t size) {
    FP *out_cope = new FP[size+1];

    cope->extend64(out_cope, size+1, FP::slot_stride_bits, num_pack);
    sender_check64(out_cope, out_cope[size], size);

    for(std::size_t i = 0; i < size; ++i)
      share[i] = out_cope[i];

    delete[] out_cope;
  }

  // recver
  template<typename FPS>
  void compute_recv64(FPS *share, FP *x, std::size_t size) {
    FP *out_cope = new FP[size+1];

    cope->extend64(out_cope, x, size+1, FP::slot_stride_bits, num_pack);
    recver_check64(out_cope, x, out_cope[size], x[size], size);

    for(std::size_t i = 0; i < size; ++i) {
      share[i].setLowHigh(out_cope[i], x[i]);
    }

    delete[] out_cope;
  }

  // recver
  template<typename FPS>
  void compute_recv64(FPS *share, std::size_t size) {
    uint64_t *buf = new uint64_t[size+1];
    prg.random_data_unaligned(buf, (size+1) * sizeof(uint64_t));
    FP *x = new FP[size+1];
    if(FP::PR_num_pack == 1) {
      for(std::size_t i = 0; i < size+1; ++i)
        x[i] = buf[i];
    } else {
      for(std::size_t i = 0; i < size+1; ++i) {
        buf[i] = buf[i] & FP::PR_mask;
        x[i].assign_no_mod(FP::copy_compose(buf[i]));
      }
    }
    delete[] buf;

    compute_recv64(share, x, size);
    delete[] x;
  }

  // sender check
  void sender_check64(FP *share, FP b, std::size_t size) {
    FP seed;
    seed.rand(prg);
    seed.send(io);
    FP *chi = new FP[size];
    check_coeff_gen(chi, seed, size);

    FP y = field_inn_prdt_sum_red(share, chi, size);
    y = y + b;
    FP xz[2];
    io->recv_data(xz, 2 * sizeof(FP));
    xz[1] = xz[1] * Delta;
    y = y + xz[1];
    if (y != xz[0]) {
      error("base sVOLE check fails");
    }
    delete[] chi;
  }

  // receiver check
  void recver_check64(FP *share, FP *x, FP c, FP a,
                    std::size_t size) {
    FP seed;
    seed.recv(io);
    FP *chi = new FP[size];
    check_coeff_gen(chi, seed, size);

    FP xz[2];
    xz[0] = field_inn_prdt_sum_red(share, chi, size);
    xz[1] = field_inn_prdt_sum_red(x, chi, size);
    xz[0] = xz[0] + c;
    xz[1] = xz[1] + a;
    io->send_data(xz, 2 * sizeof(FP));
    delete[] chi;
  }
};

#endif
