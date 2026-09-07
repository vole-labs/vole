#ifndef VOLE_MVOLE_PRIMAL_H__
#define VOLE_MVOLE_PRIMAL_H__

// Primal-LPN n-party VOLE with a final consistency check (Protocol Pi_nVOLE of
// the MVZK paper, after Le Mans). This is the *uncommitted* alternative to
// MCVoleFp (vole/mcvole.h): the same n-party correlation, built on the plain
// primal VoleTriple instead of the dual-LPN CVoleFp, so each pairwise VOLE
// runs at the primal rate. The two classes are independent; nothing here is
// shared with the committed path.
//
// Correlation: n verifiers V_0..V_{n-1} plus a prover / king P (id n). V_i
// holds a value vector u^i, and for every j != i a MAC w^i_j = u^i * Delta_j +
// v^j_i, where V_j holds the key v^j_i and its global key Delta_j. P learns
// every u^i (from the seeds) and nothing else.
//
// Why primal works without a commitment: a VoleTriple receiver's values are a
// function of its programming seed alone (base values, MPFSS positions and
// point values all come from prog_prg; the LPN matrix is public), so V_i
// seeded with the same seed towards every peer gets the same u^i in all its
// pairwise VOLEs, and P regenerates u^i locally (PrimalValueLocal below mirrors
// the receiver pipeline exactly). What a malicious V_i could do is use
// different seeds with different peers; the consistency check catches that.
//
// Consistency check (paper Protocol 1, "Consistency Check"), run by the
// verifiers only, once per check() over all extend rounds since the last one:
//   per round, right after the pairwise extends, a coin gives chi and every
//   party folds its outputs: u = sum_l chi_l u_l + u_mask, likewise the MACs
//   and keys, where u_mask is one output per round reserved as the paper's
//   u_{m+1} so that opening the fold reveals nothing about the outputs;
//   then (steps 3-7) each V_i shares a zero-sharing b^i, all open u_hat =
//   sum_i u^i masked by the b's, commit to (u^i, {Z^i_j}) with Z^i_j = w^i_j
//   for j != i and Z^i_i = (u^i - u_hat) Delta_i - sum_{j != i} v^i_j, agree
//   on a hash of the broadcast messages, open, and check
//     u_hat == sum_i u^i,   Z^j_k == v^k_j + u^j Delta_k (all j != k),
//     sum_i Z^i_j == 0 (all j).
// The last equation is Delta_j * (sum_{i != j} u^i_{as seen by j} + u^j -
// u_hat) and is what pins every verifier's values to a single vector.
//
// Network: every verifier holds an IO** to every other party (ios[j]); the
// pairwise VOLEs use MCVoleFp's antipode-ring schedule (deadlock-free
// all-to-all: forward arc on the main thread, backward arc on a pool thread).
// The check's small broadcasts go over ios[j][0] after the mesh has joined.

#include "vole/vole-fp.h"
#include <vector>
#include <cstring>
#include <cstdlib>

using namespace emp;

// Regenerates a VoleTriple *receiver's* value vector from its programming
// seed, with no network: consumes prog_prg in exactly the order
// VoleTriple::setup / extend_recv do (choice bits before base values, per
// stage), places the MPFSS point values, and runs the same LpnFp on the value
// lane. Used by the king. Any change to VoleTriple's receiver path must be
// mirrored here (test_mvole_fp compares the two by hash).
template <typename FP, typename FPS>
class PrimalValueLocal {
public:
  PrimalLPNParameterFp61 param;
  ThreadPool *pool;
  PRG prog_prg;
  std::size_t M, ot_limit;
  std::vector<FP> pre_yz;  // base pairs for the next round (value lane)
  LpnFp<10> *lpn = nullptr;

  static constexpr bool unit_point = vole_traits<FPS>::unit_point_value;
  static constexpr std::size_t point_pairs_per_tree = unit_point ? 0 : 1;
  static constexpr std::size_t mask_pairs = vole_traits<FPS>::mask_pairs;
  static std::size_t base_pairs_for(std::size_t t) { return t * point_pairs_per_tree + mask_pairs; }

  PrimalValueLocal(ThreadPool *pool, PrimalLPNParameterFp61 param = fp_default)
      : param(param), pool(pool) {
    M = base_pairs_for(param.t) + param.k;
    ot_limit = param.n - M;
    lpn = new LpnFp<10>(param.n, param.k, pool, pool->size());
  }
  ~PrimalValueLocal() { delete lpn; }

  // Regular sparse vector: tree i gets its point at the position encoded by
  // its log_bin choice bits (SpfssRecverFp::get_index: bit 0 -> 1), value
  // base[i] (or the constant 1 for the bit instantiation).
  void place_points(FP *out, std::size_t t, std::size_t log_bin, const bool *bits, const FP *base) {
    std::size_t leave = 1ull << log_bin;
    for (std::size_t i = 0; i < t; ++i) {
      std::size_t pos = 0;
      for (std::size_t l = 0; l < log_bin; ++l) {
        pos <<= 1;
        if (!bits[i * log_bin + l]) pos += 1;
      }
      out[i * leave + pos] = unit_point ? FP(1) : base[i];
    }
  }

  void setup(block seed) {
    prog_prg.reseed(&seed);

    // ---- stage 0: choice bits, then base values (+ the base-sVOLE check mask)
    std::size_t bp0 = base_pairs_for(param.t_pre0);
    std::size_t tn0 = bp0 + param.k_pre0;
    std::size_t nb0 = param.t_pre0 * param.log_bin_sz_pre0;
    bool *bits0 = new bool[nb0];
    prog_prg.random_bool(bits0, nb0);
    std::vector<FP> x(tn0 + 1);
    for (std::size_t i = 0; i < tn0; ++i) vole_traits<FPS>::rand_base_value(x[i], prog_prg);
    x[tn0].rand(prog_prg);
    std::vector<FP> y0(param.n_pre0);
    for (auto &e : y0) e.setZero();
    place_points(y0.data(), param.t_pre0, param.log_bin_sz_pre0, bits0, x.data());
    delete[] bits0;
    LpnFp<10> lpn0(param.n_pre0, param.k_pre0, pool, pool->size());
    lpn0.compute(y0.data(), x.data() + bp0);

    // ---- stage pre: choice bits, points from stage-0 outputs
    std::size_t bp1 = base_pairs_for(param.t_pre);
    std::size_t nb1 = param.t_pre * param.log_bin_sz_pre;
    bool *bits1 = new bool[nb1];
    prog_prg.random_bool(bits1, nb1);
    std::vector<FP> y1(param.n_pre);
    for (auto &e : y1) e.setZero();
    place_points(y1.data(), param.t_pre, param.log_bin_sz_pre, bits1, y0.data());
    delete[] bits1;
    LpnFp<10> lpn1(param.n_pre, param.k_pre, pool, pool->size());
    lpn1.compute(y1.data(), y0.data() + bp1);
    pre_yz.swap(y1);
  }

  // One extend round: out has param.n entries; out[0, ot_limit) are the
  // outputs, out[ot_limit, ot_limit + M) self-seed the next round.
  void extend(FP *out) {
    std::size_t nb = param.t * param.log_bin_sz;
    bool *bits = new bool[nb];
    prog_prg.random_bool(bits, nb);
    for (std::size_t i = 0; i < param.n; ++i) out[i].setZero();
    place_points(out, param.t, param.log_bin_sz, bits, pre_yz.data());
    delete[] bits;
    lpn->compute(out, pre_yz.data() + base_pairs_for(param.t));
    pre_yz.assign(out + ot_limit, out + ot_limit + M);
  }
};

template <typename IO, typename FP, typename FPS>
class MVoleFp {
public:
  int id_party;   // 0..n_party-1 = verifiers, n_party = king / prover
  int n_party;
  std::size_t threads;
  PrimalLPNParameterFp61 param;

  FP delta;
  PRG prg;
  block prog_seed;                 // verifier: its seed
  std::vector<block> prog_seeds;   // king: every verifier's seed

  std::vector<IO **> ios;
  std::vector<VoleTriple<IO, FP, FPS> *> voleForward;  // receiver role: my u, MAC under peer's Delta
  std::vector<VoleTriple<IO, FP, FPS> *> voleInverse;  // sender role: keys for the peer's u
  std::vector<PrimalValueLocal<FP, FPS> *> locals;     // king: one per verifier
  ThreadPool *pool = nullptr;                          // backward arc of the ring
  ThreadPool *king_pool = nullptr;

  std::size_t ot_limit = 0, M = 0;
  std::size_t usable = 0;   // outputs per round handed to the caller (ot_limit - 1)
  std::size_t mask_idx = 0; // the reserved output folded as u_{m+1}

  std::vector<FPS *> fwd_buf;  // per peer, param.n entries
  std::vector<FP *> inv_buf;

  // Folded accumulators for the check (verifiers).
  FP acc_u;
  std::vector<FP> acc_mac, acc_key;
  std::size_t rounds_pending = 0;

  // Test hook: seed override for the forward VOLE towards one peer (a
  // misbehaving verifier); -1 = honest.
  int cheat_peer = -1;
  block cheat_seed;

  double t_mesh = 0, t_fold = 0, t_check = 0;  // microseconds, cumulative

  MVoleFp(int id_party, int n_party, std::size_t threads, std::vector<IO **> &ios_,
          PrimalLPNParameterFp61 param = fp_default)
      : id_party(id_party), n_party(n_party), threads(threads), param(param) {
    ios.assign(ios_.begin(), ios_.end());
    pool = new ThreadPool(1);
    if (id_party != n_party) {
      delta.rand(prg);
      if (FP::PR_num_pack > 1) delta.setHigh(FP(0));
      restrict_delta(delta);
    }
    acc_u.setZero();
    acc_mac.assign(n_party, FP());
    acc_key.assign(n_party, FP());
    for (auto &e : acc_mac) e.setZero();
    for (auto &e : acc_key) e.setZero();
  }

  ~MVoleFp() {
    for (auto p : voleForward) if (p) delete p;
    for (auto p : voleInverse) if (p) delete p;
    for (auto p : locals) if (p) delete p;
    for (auto p : fwd_buf) if (p) delete[] p;
    for (auto p : inv_buf) if (p) delete[] p;
    if (pool) delete pool;
    if (king_pool) delete king_pool;
  }

  bool is_king() const { return id_party == n_party; }

  static void fail(const char *what) {
    fprintf(stderr, "mvole: %s\n", what);
    fflush(stderr);
    std::exit(1);
  }

  // ------------------------------------------------------------------ setup
  void setup() {
    if (is_king()) {
      prog_seeds.resize(n_party);
      prg.random_block(prog_seeds.data(), n_party);
      for (int i = 0; i < n_party; ++i) {
        ios[i][0]->send_data(&prog_seeds[i], sizeof(block));
        ios[i][0]->flush();
      }
      king_pool = new ThreadPool(threads);
      locals.resize(n_party, nullptr);
      for (int i = 0; i < n_party; ++i) {
        locals[i] = new PrimalValueLocal<FP, FPS>(king_pool, param);
        locals[i]->setup(prog_seeds[i]);
      }
      ot_limit = locals[0]->ot_limit;
      M = locals[0]->M;
    } else {
      ios[n_party][0]->recv_data(&prog_seed, sizeof(block));
      voleForward.assign(n_party, nullptr);
      voleInverse.assign(n_party, nullptr);
      fwd_buf.assign(n_party, nullptr);
      inv_buf.assign(n_party, nullptr);
      mesh_run(/*is_setup=*/true);
      int j0 = (id_party == 0) ? 1 : 0;
      ot_limit = voleForward[j0]->ot_limit;
      M = voleForward[j0]->M;
      for (int j = 0; j < n_party; ++j) {
        if (j == id_party) continue;
        fwd_buf[j] = new FPS[param.n];
        inv_buf[j] = new FP[param.n];
      }
    }
    usable = ot_limit - 1;
    mask_idx = ot_limit - 1;
  }

  // ---- antipode ring schedule (same as MCVoleFp) ----
  void mesh_run(bool is_setup) {
    int half = n_party / 2;
    int dest_back = (id_party >= half) ? (id_party - half) : (id_party + half);
    int dest_frnt = (dest_back == 0) ? (n_party - 1) : (dest_back - 1);

    std::vector<std::future<void>> fut;
    fut.push_back(pool->enqueue([this, dest_back, is_setup]() {
      int stop = (dest_back == 0) ? (n_party - 1) : (dest_back - 1);
      int j = (id_party == 0) ? (n_party - 1) : (id_party - 1);
      while (j != stop) {
        if (is_setup) vole_setup_with(j); else vole_extend_with(j);
        j = (j == 0) ? (n_party - 1) : (j - 1);
      }
    }));
    int stop = (dest_frnt == n_party - 1) ? 0 : (dest_frnt + 1);
    int j = (id_party == n_party - 1) ? 0 : (id_party + 1);
    while (j != stop) {
      if (is_setup) vole_setup_with(j); else vole_extend_with(j);
      j = (j == n_party - 1) ? 0 : (j + 1);
    }
    for (auto &f : fut) f.get();
  }

  // Both sub-VOLEs with peer j; the higher id runs its receiver role first so
  // the two sides agree on the order.
  void vole_setup_with(int j) {
    if (id_party > j) { make_forward(j); make_inverse(j); }
    else              { make_inverse(j); make_forward(j); }
  }
  void make_forward(int j) {
    voleForward[j] = new VoleTriple<IO, FP, FPS>(BOB, threads, ios[j], param);
    voleForward[j]->setup_prog(j == cheat_peer ? cheat_seed : prog_seed);
    voleForward[j]->template setup<FP>();
    ios[j][0]->flush();
  }
  void make_inverse(int j) {
    voleInverse[j] = new VoleTriple<IO, FP, FPS>(ALICE, threads, ios[j], param);
    voleInverse[j]->template setup<FP>(delta);
    ios[j][0]->flush();
  }
  void vole_extend_with(int j) {
    if (id_party > j) {
      voleForward[j]->extend_recv(fwd_buf[j]);
      voleInverse[j]->extend_send(inv_buf[j]);
    } else {
      voleInverse[j]->extend_send(inv_buf[j]);
      voleForward[j]->extend_recv(fwd_buf[j]);
    }
    ios[j][0]->flush();
  }

  // ----------------------------------------------------------------- extend
  // King: one round of every verifier's values; u[i] has >= usable entries.
  void extend_local(std::vector<FP *> &u) {
    std::vector<FP> tmp(param.n);
    for (int i = 0; i < n_party; ++i) {
      locals[i]->extend(tmp.data());
      std::copy(tmp.begin(), tmp.begin() + usable, u[i]);
    }
  }

  // Verifier: one round. Outputs `usable` entries each: u (my values), mac[j]
  // (MAC on u under Delta_j), key[j] (key for peer j's values). Entries for
  // j == id_party are ignored. Then folds the round into the check
  // accumulators (one coin flip among the verifiers).
  void extend(FP *u, std::vector<FP *> &mac, std::vector<FP *> &key) {
    auto t0 = clock_start();
    mesh_run(/*is_setup=*/false);
    t_mesh += time_from(t0);
    t0 = clock_start();

    int j0 = (id_party == 0) ? 1 : 0;
    for (std::size_t l = 0; l < usable; ++l) u[l] = fwd_buf[j0][l].getHigh();
    for (int j = 0; j < n_party; ++j) {
      if (j == id_party) continue;
      for (std::size_t l = 0; l < usable; ++l) mac[j][l] = fwd_buf[j][l].getLow();
      std::copy(inv_buf[j], inv_buf[j] + usable, key[j]);
    }
    fold_round(u, mac, key);
    t_fold += time_from(t0);
  }

  // chi from a fresh coin (after the outputs are fixed), fold m outputs plus
  // the reserved mask output.
  void fold_round(const FP *u, std::vector<FP *> &mac, std::vector<FP *> &key) {
    block seed = coin();
    FP d; d.from_block(seed);
    FP *chi = new FP[usable];
    check_coeff_gen(chi, d, (int)usable);
    int j0 = (id_party == 0) ? 1 : 0;
    acc_u = acc_u + field_inn_prdt_sum_red(chi, u, usable) + fwd_buf[j0][mask_idx].getHigh();
    for (int j = 0; j < n_party; ++j) {
      if (j == id_party) continue;
      acc_mac[j] = acc_mac[j] + field_inn_prdt_sum_red(chi, mac[j], usable) + fwd_buf[j][mask_idx].getLow();
      acc_key[j] = acc_key[j] + field_inn_prdt_sum_red(chi, key[j], usable) + inv_buf[j][mask_idx];
    }
    delete[] chi;
    ++rounds_pending;
  }

  // ------------------------------------------------------------------ check
  // Protocol 1 steps 3-7 on the accumulated folds. Aborts the process on any
  // failure (every honest verifier detects it in the same step).
  void check() {
    if (is_king() || rounds_pending == 0) return;
    auto t0 = clock_start();
    const int n = n_party;

    // 3. zero-sharing b^i: b_j to V_j, own share = -sum.
    std::vector<FP> b(n), b_recv(n);
    FP sum; sum.setZero();
    for (int j = 0; j < n; ++j) {
      if (j == id_party) continue;
      b[j].rand(prg);
      sum = sum + b[j];
    }
    b[id_party] = sum.negate();
    for (int j = 0; j < n; ++j) {
      if (j == id_party) continue;
      ios[j][0]->send_data(&b[j], sizeof(FP));
      ios[j][0]->flush();
    }
    for (int j = 0; j < n; ++j) {
      if (j == id_party) continue;
      ios[j][0]->recv_data(&b_recv[j], sizeof(FP));
    }

    // 4. shares of u_hat: u^i + sum_k b^k_i, broadcast, reconstruct.
    std::vector<FP> share(n);
    share[id_party] = acc_u + b[id_party];
    for (int j = 0; j < n; ++j) if (j != id_party) share[id_party] = share[id_party] + b_recv[j];
    broadcast(share.data());
    FP u_hat; u_hat.setZero();
    for (int j = 0; j < n; ++j) u_hat = u_hat + share[j];

    // 5. commit to (u^i, Z^i_0..Z^i_{n-1}).
    std::vector<FP> Z(n);
    FP zi = (acc_u - u_hat) * delta;
    for (int j = 0; j < n; ++j) {
      if (j == id_party) continue;
      Z[j] = acc_mac[j];
      zi = zi - acc_key[j];
    }
    Z[id_party] = zi;
    std::vector<block> nonce(n), com(n);
    prg.random_block(&nonce[id_party], 1);
    com[id_party] = commit(nonce[id_party], acc_u, Z.data());
    broadcast(com.data());

    // 6. everyone agrees on the broadcast transcript (shares + commitments).
    Hash h;
    h.put(share.data(), n * sizeof(FP));
    h.put_block(com.data(), n);
    block dig = digest_block(h);
    std::vector<block> digs(n);
    digs[id_party] = dig;
    broadcast(digs.data());
    for (int j = 0; j < n; ++j)
      if (!cmpBlock(&digs[j], &dig, 1)) fail("consistency check: transcript hashes differ");

    // 7. open and verify.
    std::vector<FP> u_all(n);
    std::vector<std::vector<FP>> Z_all(n, std::vector<FP>(n));
    u_all[id_party] = acc_u;
    Z_all[id_party] = Z;
    std::vector<char> pkt(sizeof(block) + (n + 1) * sizeof(FP));
    std::memcpy(pkt.data(), &nonce[id_party], sizeof(block));
    std::memcpy(pkt.data() + sizeof(block), &acc_u, sizeof(FP));
    std::memcpy(pkt.data() + sizeof(block) + sizeof(FP), Z.data(), n * sizeof(FP));
    for (int j = 0; j < n; ++j) {
      if (j == id_party) continue;
      ios[j][0]->send_data(pkt.data(), pkt.size());
      ios[j][0]->flush();
    }
    for (int j = 0; j < n; ++j) {
      if (j == id_party) continue;
      ios[j][0]->recv_data(pkt.data(), pkt.size());
      std::memcpy(&nonce[j], pkt.data(), sizeof(block));
      std::memcpy(&u_all[j], pkt.data() + sizeof(block), sizeof(FP));
      std::memcpy(Z_all[j].data(), pkt.data() + sizeof(block) + sizeof(FP), n * sizeof(FP));
      block c = commit(nonce[j], u_all[j], Z_all[j].data());
      if (!cmpBlock(&c, &com[j], 1)) fail("consistency check: commitment opening mismatch");
    }
    FP s; s.setZero();
    for (int j = 0; j < n; ++j) s = s + u_all[j];
    if (!(s == u_hat)) fail("consistency check: u_hat != sum of opened u^i");
    for (int j = 0; j < n; ++j) {
      if (j == id_party) continue;
      FP expect = acc_key[j] + u_all[j] * delta;
      if (!(Z_all[j][id_party] == expect)) fail("consistency check: MAC/key relation fails");
    }
    for (int j = 0; j < n; ++j) {
      FP col; col.setZero();
      for (int i = 0; i < n; ++i) col = col + Z_all[i][j];
      FP zero; zero.setZero();
      if (!(col == zero)) fail("consistency check: column sum of Z is nonzero");
    }

    acc_u.setZero();
    for (auto &e : acc_mac) e.setZero();
    for (auto &e : acc_key) e.setZero();
    rounds_pending = 0;
    t_check += time_from(t0);
  }

  // ----------------------------------------------------------------- helpers
  // Broadcast slot v[id_party] to every verifier and fill the other slots.
  template <typename T>
  void broadcast(T *v) {
    for (int j = 0; j < n_party; ++j) {
      if (j == id_party) continue;
      ios[j][0]->send_data(&v[id_party], sizeof(T));
      ios[j][0]->flush();
    }
    for (int j = 0; j < n_party; ++j) {
      if (j == id_party) continue;
      ios[j][0]->recv_data(&v[j], sizeof(T));
    }
  }

  // n-party coin among the verifiers: commit to a random block, open, XOR.
  block coin() {
    std::vector<block> r(n_party), c(n_party);
    prg.random_block(&r[id_party], 1);
    c[id_party] = Hash::hash_for_block(&r[id_party], sizeof(block));
    broadcast(c.data());
    broadcast(r.data());
    block out = zero_block;
    for (int j = 0; j < n_party; ++j) {
      if (j != id_party) {
        block cj = Hash::hash_for_block(&r[j], sizeof(block));
        if (!cmpBlock(&cj, &c[j], 1)) fail("coin: opening mismatch");
      }
      out = out ^ r[j];
    }
    return out;
  }

  // SHA-256 digest truncated to one block (Hash::digest writes 32 bytes).
  static block digest_block(Hash &h) {
    alignas(block) char d[Hash::DIGEST_SIZE];
    h.digest(d);
    block out;
    std::memcpy(&out, d, sizeof(block));
    return out;
  }

  block commit(block nonce, const FP &u, const FP *Z) {
    Hash h;
    h.put_block(&nonce, 1);
    h.put(&u, sizeof(FP));
    h.put(Z, n_party * sizeof(FP));
    return digest_block(h);
  }
};

#endif  // VOLE_MVOLE_PRIMAL_H__
