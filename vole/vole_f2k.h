#ifndef VOLE_F2K_FERRET_STYLE_H__
#define VOLE_F2K_FERRET_STYLE_H__

// Ferret-COT-style interface on top of the primal silent VOLE (VoleTriple),
// modelled on emp-ot 0.3.0's FerretCOT::rcot / send_cot / recv_cot:
//
//   * one round-sized internal buffer with a cursor; any request length is
//     served from the buffer and the buffer is refilled by extending a round
//     (so callers never size buffers or count rounds themselves);
//   * random VOLE:  rvole_send(K, n) / rvole_recv(vm, n)  give K and (u, M)
//     with M = K + u * Delta, u random;
//   * chosen-input VOLE (the analogue of send_cot / recv_cot, paper App. A.9):
//     the receiver derandomises its random u to a chosen x with one message
//     f = u - x, and the sender corrects its key K' = K + f * Delta, so that
//     M = K' + x * Delta.
//
// Roles follow VoleTriple: ALICE holds Delta and the keys ("sender", as in
// Ferret's COT sender); BOB holds the values and MACs.
//
// RVole<IO, FP, FPS> is field-generic (any of the vole/fields types).
// F2kVole<IO> specialises it to GF(2^128) and exposes plain `block` arguments,
// the way FerretCOT does, so it drops into code written against emp-ot's F2k
// VOLE / Ferret shape without touching the field classes.

#include "vole/vole-fp.h"
#include "vole/fields/fp2x128.h"
#include "vole/fields/f2.h"
#include <cstring>
#include <vector>

template <typename IO, typename FP, typename FPS>
class RVole {
public:
  int party;
  VoleTriple<IO, FP, FPS> vt;
  IO *io;
  FP Delta;

  // Internal round buffer (param.n entries; the first ot_limit are usable,
  // the tail is the self-seeding reserve that VoleTriple copies out itself).
  FP *buf_key = nullptr;    // ALICE
  FPS *buf_vm = nullptr;    // BOB
  std::size_t used = 0;     // cursor into the usable prefix
  bool ready = false;

  RVole(int party, std::size_t threads, IO **ios,
        PrimalLPNParameterFp61 param = fp_default)
      : party(party), vt(party, threads, ios, param), io(ios[0]) {}

  ~RVole() {
    if (buf_key) delete[] buf_key;
    if (buf_vm) delete[] buf_vm;
  }

  std::size_t ot_limit() const { return vt.ot_limit; }
  std::size_t left() const { return vt.ot_limit - used; }

  // ---- setup ----
  void setup(FP delta) {  // ALICE
    Delta = delta;
    vt.template setup<FP>(delta);
    buf_key = new FP[vt.param.n];
    used = vt.ot_limit;  // empty
    ready = true;
  }
  void setup() {          // BOB
    PRG prg;
    block s;
    prg.random_block(&s, 1);
    vt.setup_prog(s);
    vt.template setup<FP>();
    buf_vm = new FPS[vt.param.n];
    used = vt.ot_limit;
    ready = true;
  }

  // ---- random VOLE, any length (Ferret's rcot) ----
  void rvole_send(FP *key, std::size_t num) {
    if (!ready || party != ALICE) error("RVole::rvole_send: not set up as ALICE");
    while (num > 0) {
      if (left() == 0) { vt.extend_send(buf_key); used = 0; }
      std::size_t take = std::min(num, left());
      std::copy(buf_key + used, buf_key + used + take, key);
      key += take; used += take; num -= take;
    }
  }
  void rvole_recv(FPS *vm, std::size_t num) {
    if (!ready || party != BOB) error("RVole::rvole_recv: not set up as BOB");
    while (num > 0) {
      if (left() == 0) { vt.extend_recv(buf_vm); used = 0; }
      std::size_t take = std::min(num, left());
      std::copy(buf_vm + used, buf_vm + used + take, vm);
      vm += take; used += take; num -= take;
    }
  }

  // ---- chosen-input VOLE (Ferret's send_cot / recv_cot) ----
  // BOB chooses x; afterwards M = K' + x * Delta.
  void vole_send(FP *key, std::size_t num) {  // ALICE
    rvole_send(key, num);
    std::vector<FP> f(num);
    io->recv_data(f.data(), num * sizeof(FP));
    for (std::size_t i = 0; i < num; ++i)
      key[i] = key[i] + FP(f[i].val) * Delta;  // canonical f
  }
  void vole_recv(const FP *x, FP *mac, std::size_t num) {  // BOB
    std::vector<FPS> vm(num);
    rvole_recv(vm.data(), num);
    std::vector<FP> f(num);
    for (std::size_t i = 0; i < num; ++i) {
      f[i] = vm[i].getHigh() - x[i];  // u - x
      mac[i] = vm[i].getLow();
    }
    io->send_data(f.data(), num * sizeof(FP));
    io->flush();
  }
};

// GF(2^128) instance with block-typed arguments (Ferret / emp-ot F2kVOLE shape).
template <typename IO>
class F2kVole : public RVole<IO, FP2x128, FP2x128x2> {
public:
  using Base = RVole<IO, FP2x128, FP2x128x2>;
  using Base::party;

  F2kVole(int party, std::size_t threads, IO **ios,
          PrimalLPNParameterFp61 param = fp_default)
      : Base(party, threads, ios, param) {}

  static block to_block(const FP2x128 &v) { return FP2x128::to_block(v.val); }
  static FP2x128 from_block(block b) { FP2x128 r; r.val = FP2x128::from_blk(b); return r; }

  // ALICE: random Delta unless given (like FerretCOT, which samples its own).
  void setup(block delta) { Base::setup(from_block(delta)); }
  void setup() {
    if (party == ALICE) {
      PRG prg; block d; prg.random_block(&d, 1);
      Base::setup(from_block(d));
    } else {
      Base::setup();
    }
  }
  block delta() const { return to_block(this->Delta); }

  // Random VOLE: ALICE gets keys; BOB gets (val, mac) with mac = key ^ gfmul(val, Delta).
  void rvole(block *key, std::size_t num) {                 // ALICE
    std::vector<FP2x128> k(num);
    Base::rvole_send(k.data(), num);
    for (std::size_t i = 0; i < num; ++i) key[i] = to_block(k[i]);
  }
  void rvole(block *val, block *mac, std::size_t num) {     // BOB
    std::vector<FP2x128x2> vm(num);
    Base::rvole_recv(vm.data(), num);
    for (std::size_t i = 0; i < num; ++i) {
      val[i] = FP2x128::to_block(vm[i].val[1]);
      mac[i] = FP2x128::to_block(vm[i].val[0]);
    }
  }

  // Chosen-input VOLE: BOB supplies x; afterwards mac = key ^ gfmul(x, Delta).
  void vole_send(block *key, std::size_t num) {             // ALICE
    std::vector<FP2x128> k(num);
    Base::vole_send(k.data(), num);
    for (std::size_t i = 0; i < num; ++i) key[i] = to_block(k[i]);
  }
  void vole_recv(const block *x, block *mac, std::size_t num) {  // BOB
    std::vector<FP2x128> xf(num), m(num);
    for (std::size_t i = 0; i < num; ++i) xf[i] = from_block(x[i]);
    Base::vole_recv(xf.data(), m.data(), num);
    for (std::size_t i = 0; i < num; ++i) mac[i] = to_block(m[i]);
  }
};

// F_2-value instantiation (Ferret's correlated OT): receiver blocks carry the
// choice bit in the LSB, R = K xor b*Delta. Same rcot / send_cot / recv_cot
// surface as emp-ot's FerretCOT, including 1-bit derandomisation.
template <typename IO>
class FerretVole : public RVole<IO, F2kKey, F2Auth> {
public:
  using Base = RVole<IO, F2kKey, F2Auth>;
  using Base::party;
  using Base::io;

  FerretVole(int party, std::size_t threads, IO **ios,
             PrimalLPNParameterFp61 param = fp_ferret_f2)
      : Base(party, threads, ios, param) {}

  // ALICE: Delta with LSB 1 (random unless given); BOB: no argument.
  void setup(block delta) {
    F2kKey d; d.val = FP2x128::from_blk(delta); restrict_delta(d);
    Base::setup(d);
  }
  void setup() {
    if (party == ALICE) {
      PRG prg; block d; prg.random_block(&d, 1);
      setup(d);
    } else {
      Base::setup();
    }
  }
  block delta() const { return FP2x128::to_block(this->Delta.val); }

  // Random COT: ALICE keys (LSB 0), BOB blocks with LSB = choice bit.
  void rcot(block *data, std::size_t num) {
    if (party == ALICE) {
      std::vector<F2kKey> k(num);
      Base::rvole_send(k.data(), num);
      for (std::size_t i = 0; i < num; ++i) data[i] = FP2x128::to_block(k[i].val);
    } else {
      std::vector<F2Auth> r(num);
      Base::rvole_recv(r.data(), num);
      for (std::size_t i = 0; i < num; ++i) data[i] = r[i].to_block();
    }
  }

  // Chosen choice bits, Ferret's online phase: BOB sends d_i = LSB(R_i) xor
  // b_i (one bit each) and keeps R_i; ALICE sets K_i' = K_i xor d_i*Delta, so
  // R_i = K_i' xor b_i*Delta.
  void send_cot(block *data, std::size_t num) {          // ALICE
    rcot(data, num);
    bool *d = new bool[num];
    io->recv_data(d, num);
    block Dl = delta();
    for (std::size_t i = 0; i < num; ++i) if (d[i]) data[i] = data[i] ^ Dl;
    delete[] d;
  }
  void recv_cot(block *data, const bool *b, std::size_t num) {  // BOB
    rcot(data, num);
    bool *d = new bool[num];
    for (std::size_t i = 0; i < num; ++i) d[i] = getLSB(data[i]) != b[i];
    io->send_data(d, num);
    io->flush();
    delete[] d;
  }
};

#endif  // VOLE_F2K_FERRET_STYLE_H__
