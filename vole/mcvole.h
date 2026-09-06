#ifndef DIST_PSI_MCVOLE_FP_H__
#define DIST_PSI_MCVOLE_FP_H__

// Field-generic n-party *pairwise committed* VOLE (the multi-party VOLE of
// paper Sec. 6.2, as used by multi-verifier ZK [EPSW24]), built on
// CVoleFp. One extra party P_king picks a seed for every one of the n parties,
// reproduces each party's committed value/commitment locally (CVoleFp's
// king-local mode), and publishes the commitments {com^i} to everyone. Each
// ordered pair (i, j) then runs a committed VOLE: i is the committer (its seed
// fixes u^i and com^i), j is the verifier. Because com^i is public, every
// verifier checks its key share against the published commitment, so the usual
// SPDZ MAC check at the end is unnecessary.
//
// Parties: ids 0..n-1 are the n VOLE parties (verifiers/committers); id == n is
// P_king. Each party holds an IO** channel to every other party (ios[j]).
//
// Self-consistency of u^i across i's committer instances: CVoleFp always
// materializes the value and self-seeds from it, and every voleForward[j] of
// party i is seeded with the SAME prog_seed, so they all reproduce the SAME
// u^i / com^i (only the MACs differ per verifier j's Delta). No explicit value
// propagation is needed.
//
// Network topology / scheduling (antipode-ring schedule): the
// all-to-all mesh is split by the antipode ring. Party id processes a forward
// arc {id+1 .. dest_frnt} on the main thread and a backward arc {id-1 .. dest_back}
// on a pool thread, with dest_back = id +/- n/2 and dest_frnt = dest_back - 1.
// For any pair, one side runs it on its forward (main) thread and the other on
// its backward (pool) thread, so the two halves are step-synchronized and the
// mesh cannot deadlock.

#include "vole/cvole.h"

using namespace emp;

template <typename IO, typename FP, typename FPS>
class MCVoleFp {
public:
  int id_party;          // 0..n_party-1 = VOLE parties, n_party = king
  int n_party;           // number of VOLE parties (excludes the king)
  std::size_t threads;
  std::size_t target_size;
  CVoleFpParam param;

  FP delta;              // this party's VOLE key (verifier role)
  PRG prg;

  std::size_t rounds = 0;  // extend rounds to cover target_size
  std::size_t cn = 0;      // commitment buffer length per party (= rounds*n_com)

  std::vector<block> prog_seeds;               // king: a seed per VOLE party
  block prog_seed;                             // VOLE party: its own seed
  std::vector<std::vector<FP>> published_com;  // com^i for every i (length cn)

  std::vector<IO **> ios;                                // ios[j] -> party j
  std::vector<CVoleFp<IO, FP, FPS> *> voleForward;       // committer role (per peer)
  std::vector<CVoleFp<IO, FP, FPS> *> voleInverse;       // verifier role (per peer)

  ThreadPool *pool = nullptr;  // runs the backward arc of the ring

  MCVoleFp(int id_party, int n_party, std::size_t threads,
           std::vector<IO **> &ios_,
           std::size_t target_size = cvole_fp_default.n)
      : id_party(id_party), n_party(n_party), threads(threads),
        target_size(target_size) {
    param = cvole_resolve_param<FP>(cvole_fp_param_for(target_size));
    std::size_t M = (param.t + 1) + (param.t_com + 1);
    std::size_t ot_limit = param.n - M;
    rounds = (target_size + ot_limit - 1) / ot_limit;
    cn = rounds * param.n_com;
    ios.assign(ios_.begin(), ios_.end());
    pool = new ThreadPool(1);
    if (id_party != n_party) {
      delta.rand(prg);
      if (FP::PR_num_pack > 1) delta.setHigh(FP(0));
      restrict_delta(delta);  // ring: clamp Delta to its s-bit sub-ring (no-op for fields)
    }
  }

  ~MCVoleFp() {
    if (pool != nullptr) delete pool;
    for (auto p : voleForward) if (p != nullptr) delete p;
    for (auto p : voleInverse) if (p != nullptr) delete p;
  }

  bool is_king() const { return id_party == n_party; }

  // ---- King: distribute seeds, then publish each party's commitment ----
  void run_king() {
    prog_seeds.resize(n_party);
    prg.random_block(prog_seeds.data(), n_party);
    for (int i = 0; i < n_party; ++i) {
      ios[i][0]->send_data(&prog_seeds[i], sizeof(block));
      ios[i][0]->flush();
    }

    // Reproduce com^i locally from each seed (no interaction).
    published_com.resize(n_party);
    for (int i = 0; i < n_party; ++i) {
      CVoleFp<IO, FP, FPS> local(threads, param);
      local.setup_prog(prog_seeds[i]);
      local.setup_local();
      std::size_t xn = local.x_buf_size(rounds);
      FP *x_loc = new FP[xn];
      FP *com_loc = new FP[cn];
      local.extend_inplace_local(x_loc, com_loc, rounds);
      published_com[i].assign(com_loc, com_loc + cn);
      delete[] x_loc;
      delete[] com_loc;
    }

    // Publish the full commitment set to every party.
    for (int j = 0; j < n_party; ++j) {
      for (int i = 0; i < n_party; ++i)
        ios[j][0]->send_data(published_com[i].data(), cn * sizeof(FP));
      ios[j][0]->flush();
    }
  }

  // ---- VOLE party: receive seed + commitments, run the pairwise mesh ----
  void run_party() {
    ios[n_party][0]->recv_data(&prog_seed, sizeof(block));
    published_com.resize(n_party);
    for (int i = 0; i < n_party; ++i) {
      published_com[i].resize(cn);
      ios[n_party][0]->recv_data(published_com[i].data(), cn * sizeof(FP));
    }

    voleForward.resize(n_party, nullptr);
    voleInverse.resize(n_party, nullptr);

    mesh_run(/*is_setup=*/true);
    mesh_run(/*is_setup=*/false);
  }

  void run() {
    if (is_king()) run_king();
    else run_party();
  }

  // ---- antipode ring schedule (deadlock-free all-to-all) ----
  void mesh_run(bool is_setup) {
    int half = n_party / 2;
    int dest_back = (id_party >= half) ? (id_party - half) : (id_party + half);
    int dest_frnt = (dest_back == 0) ? (n_party - 1) : (dest_back - 1);

    // backward arc on the pool thread: id-1, id-2, ..., dest_back (inclusive)
    std::vector<std::future<void>> fut;
    fut.push_back(pool->enqueue([this, dest_back, is_setup]() {
      int stop = (dest_back == 0) ? (n_party - 1) : (dest_back - 1);
      int j = (id_party == 0) ? (n_party - 1) : (id_party - 1);
      while (j != stop) {
        if (is_setup) vole_setup_with(j);
        else vole_extend_with(j);
        j = (j == 0) ? (n_party - 1) : (j - 1);
      }
    }));

    // forward arc on the main thread: id+1, id+2, ..., dest_frnt (inclusive)
    int stop = (dest_frnt == n_party - 1) ? 0 : (dest_frnt + 1);
    int j = (id_party == n_party - 1) ? 0 : (id_party + 1);
    while (j != stop) {
      if (is_setup) vole_setup_with(j);
      else vole_extend_with(j);
      j = (j == n_party - 1) ? 0 : (j + 1);
    }
    for (auto &f : fut) f.get();
  }

  // Per-peer schedule trace, gated on MCVOLE_DEBUG (a mesh deadlock is most
  // easily diagnosed by seeing which peer each party last engaged).
  void dbg(const char *what, int j) {
    if (getenv("MCVOLE_DEBUG"))
      fprintf(stderr, "[p%d] %s peer %d\n", id_party, what, j), fflush(stderr);
  }

  // Construct + set up both sub-VOLEs with peer j. The higher-id party runs its
  // committer sub-VOLE first so the two parties agree on the order:
  //   VOLE-1: max-id committer  <->  min-id verifier
  //   VOLE-2: max-id verifier   <->  min-id committer
  void vole_setup_with(int j) {
    dbg("setup BEGIN", j);
    if (id_party > j) {
      make_forward(j);
      make_inverse(j);
    } else {
      make_inverse(j);
      make_forward(j);
    }
    dbg("setup END", j);
  }

  void make_forward(int j) {  // this party = committer, peer j = verifier
    voleForward[j] = new CVoleFp<IO, FP, FPS>(BOB, threads, ios[j], param);
    voleForward[j]->setup_prog(prog_seed);  // before setup: reproducible value
    voleForward[j]->setup();
    // The base-sVOLE check leaves the receiver's last message buffered; flush so
    // the peer (blocked in its sender check) isn't stuck waiting on this socket
    // while we move on to a different peer.
    ios[j][0]->flush();
  }

  void make_inverse(int j) {  // this party = verifier, peer j = committer
    voleInverse[j] = new CVoleFp<IO, FP, FPS>(ALICE, threads, ios[j], param);
    voleInverse[j]->setup(delta);
    ios[j][0]->flush();
  }

  void vole_extend_with(int j) {
    dbg("extend BEGIN", j);
    if (id_party > j) {
      commit_round(j);
      verify_round(j);
    } else {
      verify_round(j);
      commit_round(j);
    }
    dbg("extend END", j);
  }

  // Committer side: produce u^i / com^i and send Hash(M[com^i]) to verifier j.
  void commit_round(int j) {
    std::size_t xn = voleForward[j]->x_buf_size(rounds);
    FPS *x_vm = new FPS[xn];
    FPS *com_vm = new FPS[cn];
    voleForward[j]->extend_inplace_recv(x_vm, com_vm, rounds);
    voleForward[j]->consistency_send_machash(com_vm, cn);
    delete[] x_vm;
    delete[] com_vm;
  }

  // Verifier side: produce key shares and check them against the king's
  // published com^j (committer j's commitment). Errors out on any mismatch.
  void verify_round(int j) {
    std::size_t xn = voleInverse[j]->x_buf_size(rounds);
    FP *x_key = new FP[xn];
    FP *com_key = new FP[cn];
    voleInverse[j]->extend_inplace_send(x_key, com_key, rounds);
    voleInverse[j]->consistency_check_pub(com_key, published_com[j].data(), cn);
    delete[] x_key;
    delete[] com_key;
  }
};

#endif  // DIST_PSI_MCVOLE_FP_H__
