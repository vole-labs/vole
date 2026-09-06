#ifndef DIST_PSI_NCVOLE_FP_H__
#define DIST_PSI_NCVOLE_FP_H__

// Multi-client committed-VOLE driver (field-generic).
//
// One server (the
// committer) runs committed VOLE with many clients (verifiers), reusing a
// single committed input x across all of them.
//
// CVoleFp's committed input is reproducible from a seed, so the server
// regenerates the *same* regular-sparse e_u/e_r -- hence the same x and the same
// commitment com -- for every client. Only the per-client MACs differ, because
// each client contributes its own Delta. Every client checks its VOLE against
// the published commitment, which proves the server uses the same input with
// all clients (the consistency guarantee crowd PSI relies on).
//
// Roles: server = CVoleFp BOB (committer, holds x and the per-client M[x]);
//        client = CVoleFp ALICE (verifier, holds its Delta and K[x]).
// id_party: 0 = server; 1..n_client = clients.

#include "vole/cvole.h"

using namespace emp;

template <typename IO, typename FP, typename FPS>
class ProgNCVoleFp {
public:
  int id_party;
  int n_client;
  std::size_t threads;
  CVoleFpParam param;
  std::size_t target_size;  // number of committed VOLE outputs needed

  block committer_seed;  // server only: fixes the committed input
  FP delta;              // client only: this client's VOLE key
  PRG prg;

  // server only: the published commitment value, materialized at the first
  // client and reused for every later client's consistency check.
  std::vector<FP> committed_com_val;
  bool com_val_ready = false;

  std::vector<CVoleFp<IO, FP, FPS> *> cvoles;  // server: one per client
  CVoleFp<IO, FP, FPS> *cvole = nullptr;       // client

  // (#5) the per-round LPN param is selected to cover target_size; the number
  // of rounds (and buffer sizes) depends on ot_limit, known only after setup.
  ProgNCVoleFp(int id_party, int n_client, std::size_t threads,
               std::size_t target_size = cvole_fp_default.n)
      : id_party(id_party), n_client(n_client), threads(threads),
        target_size(target_size) {
    param = cvole_resolve_param<FP>(cvole_fp_param_for(target_size));
    if (id_party == 0) {
      prg.random_block(&committer_seed, 1);
      cvoles.resize(n_client, nullptr);
    }
  }

  ~ProgNCVoleFp() {
    for (auto p : cvoles)
      if (p != nullptr) delete p;
    if (cvole != nullptr) delete cvole;
  }

  // ---- server (committer) ----
  // Each client gets a fresh CVoleFp seeded with the same committer_seed, so the
  // committed input (x, com values) is reproduced identically; only the MACs
  // differ per client. The first client materializes the published commitment
  // value; later clients reuse it, and each client's own consistency check
  // (client side) enforces that it uses that one committed input.
  void extend_server(int client_id, IO **ios) {
    CVoleFp<IO, FP, FPS> *cv = new CVoleFp<IO, FP, FPS>(BOB, threads, ios, param);
    cvoles[client_id] = cv;
    cv->setup_prog(committer_seed);  // before setup: reproducible bootstrap
    cv->setup();

    std::size_t rounds = cv->rounds_for(target_size);
    std::size_t cn = cv->com_buf_size(rounds);
    FPS *x_vm = new FPS[cv->x_buf_size(rounds)];
    FPS *com_vm = new FPS[cn];
    cv->extend_inplace_recv(x_vm, com_vm, rounds);

    if (!com_val_ready) {
      committed_com_val.resize(cn);
      for (std::size_t i = 0; i < cn; ++i)
        committed_com_val[i] = com_vm[i].getHigh();
      com_val_ready = true;
    }
    cv->consistency_check_commit(committed_com_val.data(), com_vm, cn);

    delete[] x_vm;
    delete[] com_vm;
  }

  // ---- client (verifier) ----
  // Run committed VOLE with the server, obtaining keys K[x], K[com], and verify
  // K[com] + com*Delta == M[com] against the published commitment.
  void extend_client(IO **ios) {
    delta.rand(prg);
    if (FP::PR_num_pack > 1)
      delta.setHigh(FP(0));
    restrict_delta(delta);  // ring: clamp Delta to its s-bit sub-ring (no-op for fields)
    cvole = new CVoleFp<IO, FP, FPS>(ALICE, threads, ios, param);
    cvole->setup(delta);

    std::size_t rounds = cvole->rounds_for(target_size);
    std::size_t cn = cvole->com_buf_size(rounds);
    FP *x_key = new FP[cvole->x_buf_size(rounds)];
    FP *com_key = new FP[cn];
    cvole->extend_inplace_send(x_key, com_key, rounds);
    cvole->consistency_check_verify(com_key, cn);

    delete[] x_key;
    delete[] com_key;
  }
};

#endif // DIST_PSI_NCVOLE_FP_H__
