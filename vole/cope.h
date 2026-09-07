#ifndef COPE_H__
#define COPE_H__

#include "emp-ot/emp-ot.h"
#include "emp-tool/emp-tool.h"
#include "vole/fields/field_config.h"

template <typename IO, typename FP> 
class CopeFp {
public:
  int party;
  IO *io;
  block *K = nullptr;
  static constexpr std::size_t kMinBaseOT = 80;  // CSW batch minimum
  FP delta;
  PRG *G0 = nullptr, *G1 = nullptr;
  bool *delta_bool = nullptr;

  CopeFp(int party, IO *io) {
    this->party = party;
    this->io = io;
  }

  ~CopeFp() {
    if (G0 != nullptr)
      delete[] G0;
    if (G1 != nullptr)
      delete[] G1;
    if (delta_bool != nullptr)
      delete[] delta_bool;
  }

  // sender
  void initialize(FP delta, std::size_t m) {
    this->delta = delta;
    delta_bool = new bool[m];
    delta64_to_bool(delta_bool, delta.val, m);

    // CSW's extraction argument needs a batch of at least 80 OTs (its header:
    // "length must be >= 80, ~2 sigma for sigma = 40"); fp61 (61 bits) and
    // z2k (64) are below that, so pad the batch with random choice bits and
    // ignore the extra outputs on both sides.
    std::size_t m_ot = std::max<std::size_t>(m, kMinBaseOT);
    bool *choice = new bool[m_ot];
    for (std::size_t i = 0; i < m; ++i) choice[i] = delta_bool[i];
    if (m_ot > m) { PRG pad; pad.random_bool(choice + m, m_ot - m); }
    K = new block[m_ot];
    emp::CSW otco(io);  // base OT (malicious-secure); 0.3.0 used OTCO
    otco.recv(K, choice, m_ot);
    delete[] choice;

    G0 = new PRG[m];
    for (std::size_t i = 0; i < m; ++i)
      G0[i].reseed(K + i);

    delete[] K;
  }

  // recver
  void initialize(std::size_t m) {
    std::size_t m_ot = std::max<std::size_t>(m, kMinBaseOT);  // see initialize(delta, m)
    K = new block[2 * m_ot];
    PRG prg;
    prg.random_block(K, 2 * m_ot);
    emp::CSW otco(io);
    otco.send(K, K + m_ot, m_ot);

    G0 = new PRG[m];
    G1 = new PRG[m];
    for (std::size_t i = 0; i < m; ++i) {
      G0[i].reseed(K + i);
      G1[i].reseed(K + m_ot + i);
    }

    delete[] K;
  }

  // The COPE batch uses raw integer buffers as byte-pack carriers for OT
  // outputs and challenges. The buffer width must match the FIELD's storage
  // type — uint64 for FP61 (PR < 2^64) but u128 for FP107 (PR fits in 107
  // bits). Otherwise challenges and PRG outputs silently truncate to 64 bits
  // and the COPE consistency check (`b + Delta * a == c`) fails.
  using buf_t = decltype(FP().val);

  // sender batch
  void extend64(FP *ret, std::size_t size, std::size_t m, std::size_t n_pack) {
    FP *w = new FP[m * size];
    FP *v = new FP[m * size];
    buf_t *buf = new buf_t[size];
    for (std::size_t i = 0; i < m; ++i) {
      for(std::size_t j = 0; j < size; ++j) w[i*size+j].assign_no_mod((buf_t)0);
      for(std::size_t k = 0; k < n_pack; ++k) {
        G0[k*m+i].random_data_unaligned(buf, size * sizeof(buf_t));
        for (std::size_t j = 0; j < size; ++j) {
          w[i * size + j] = w[i * size + j] + (((buf[j] >> k*m)&FP::PR_mask) << (k*m));
        }
      }
    }

    FP ch[2];
    ch[0].assign_no_mod((buf_t)0);
    for (std::size_t i = 0; i < m; ++i) {
      io->recv_data(buf, size * sizeof(buf_t));
      for (std::size_t j = 0; j < size; ++j) {
        ch[1].assign_no_mod(buf[j]);
        buf_t tmp = (buf_t)0;
        for(std::size_t k = 0; k < n_pack; ++k) {
          buf_t choice = ch[delta_bool[k*m+i]].value();
          tmp |= ((choice >> (k*m)) & FP::PR_mask) << (k*m);
        }
        v[i * size + j] = w[i * size + j] + tmp;
      }
    }

    prm2pr64(ret, v, size, m);

    delete[] w;
    delete[] v;
    delete[] buf;
  }

  // recver batch
  void extend64(FP *ret, FP *u, std::size_t size, std::size_t m, std::size_t n_pack) {
    FP *w0 = new FP[m * size];
    FP *w1 = new FP[m * size];
    buf_t *buf0 = new buf_t[size];
    buf_t *buf1 = new buf_t[size];
    for (std::size_t i = 0; i < m; ++i) {
      for(std::size_t j = 0; j < size; ++j) {
        w0[i*size+j].assign_no_mod((buf_t)0);
        w1[i*size+j].assign_no_mod((buf_t)0);
      }
      for(std::size_t k = 0; k < n_pack; ++k) {
        G0[k*m+i].random_data_unaligned(buf0, size * sizeof(buf_t));
        G1[k*m+i].random_data_unaligned(buf1, size * sizeof(buf_t));
        for (std::size_t j = 0; j < size; ++j) {
          w0[i * size + j] = w0[i * size + j] + (((buf0[j] >> (k*m)) & FP::PR_mask) << (k*m));
          w1[i * size + j] = w1[i * size + j] + (((buf1[j] >> (k*m)) & FP::PR_mask) << (k*m));
        }
      }
    }

    for (std::size_t i= 0; i < m; ++i) {
      for(std::size_t j = 0; j < size; ++j) {
        w1[i * size + j] = w1[i * size + j] + u[j];
        w1[i * size + j] = w0[i * size + j] - w1[i * size + j];
        buf0[j] =  w1[i * size + j].value();
      }
      io->send_data(buf0, size * sizeof(buf_t));
    }

    prm2pr64(ret, w0, size, m);

    delete[] w0;
    delete[] w1;
    delete[] buf0;
    delete[] buf1;
  }

  // Templated so the caller's `in` keeps its full width (u128 for FP107,
  // uint64 for FP61). Otherwise `in` truncates to 64 bits and any delta bit
  // at position >=64 is silently zero, breaking the COPE correlation.
  template<typename T>
  void delta64_to_bool(bool *bdata, T in, std::size_t m) {
    for (std::size_t i = 0; i < m; ++i) {
      bdata[i] = ((in & (T)1) == (T)1);
      in >>= 1;
    }
  }

  void prm2pr64(FP *ret, FP *a, std::size_t size, std::size_t m) {
    for(std::size_t i = 0; i < size; ++i) ret[i].assign_no_mod(0ULL);
    for (std::size_t i = 0; i < m; ++i) {
      for (std::size_t j = 0; j < size; ++j) {
        FP tmp = (a[i * size + j] << i);
        ret[j] = ret[j] + tmp;
      }
    }
  }

  // debug function
  void check_triple(FP *a, FP *b, std::size_t sz) {
    if (party == ALICE) {
      io->send_data(a, sizeof(uint64_t));
      io->send_data(b, sz * sizeof(uint64_t));
    } else {
      FP delta;
      FP *c = new FP[sz];
      io->recv_data(&delta, sizeof(uint64_t));
      // Cast value() to uint64_t for printing — works for FP61 (where value
      // is already 64-bit) and surfaces only the low 64 bits for FP107.
      std::cout << "delta: " << (uint64_t)delta.value() << std::endl;
      io->recv_data(c, sz * sizeof(uint64_t));
      for (std::size_t i = 0; i < sz; ++i) {
        FP tmp = a[i] * delta;
        tmp = tmp + c[i];
        if (tmp != b[i]) {
          std::cout << "wrong triple" << i << std::endl;
          std::cout << (uint64_t)tmp.value() << " " << (uint64_t)b[i].value() << std::endl;
          abort();
        }
      }
      delete[] c;
    }
    std::cout << "pass check" << std::endl;
  }
};

#endif
