#ifndef _VOLE_TRIPLE_H_
#define _VOLE_TRIPLE_H_
#include "vole/base_cot.h"
#include "vole/base_svole.h"
#include "vole/lpn.h"
#include "vole/mpfss_reg.h"
#include <algorithm>

class PrimalLPNParameterFp61 {
public:
  std::size_t n, t, k, log_bin_sz;
  std::size_t n_pre, t_pre, k_pre, log_bin_sz_pre;
  std::size_t n_pre0, t_pre0, k_pre0, log_bin_sz_pre0;

  constexpr PrimalLPNParameterFp61() : n(0), t(0), k(0), log_bin_sz(0), n_pre(0), t_pre(0), k_pre(0),
      log_bin_sz_pre(0), n_pre0(0), t_pre0(0), k_pre0(0), log_bin_sz_pre0(0) {}
  constexpr PrimalLPNParameterFp61(std::size_t n, std::size_t t, std::size_t k, std::size_t log_bin_sz,
                         std::size_t n_pre, std::size_t t_pre, std::size_t k_pre,
                         std::size_t log_bin_sz_pre, std::size_t n_pre0, std::size_t t_pre0,
                         std::size_t k_pre0, std::size_t log_bin_sz_pre0)
      : n(n), t(t), k(k), log_bin_sz(log_bin_sz), n_pre(n_pre), t_pre(t_pre),
        k_pre(k_pre), log_bin_sz_pre(log_bin_sz_pre), n_pre0(n_pre0),
        t_pre0(t_pre0), k_pre0(k_pre0), log_bin_sz_pre0(log_bin_sz_pre0) {

    if (n != t * (1 << log_bin_sz) || n_pre != t_pre * (1 << log_bin_sz_pre) ||
        n_pre0 != t_pre0 * (1 << log_bin_sz_pre0) ||
        n_pre < k + t + 1 || n_pre0 < k_pre + t_pre + 1)
      error("LPN parameter not matched");
    // LpnFp draws idx = r & (k - 1): uniform only for a power-of-two k.
    if (!is_pow2(k) || !is_pow2(k_pre) || !is_pow2(k_pre0))
      error("LPN dimension k must be a power of two");
  }
  static constexpr bool is_pow2(std::size_t x) { return x != 0 && (x & (x - 1)) == 0; }
  std::size_t buf_sz() const { return n - t - k - 1; }
};

// Primal-LPN parameter sets: emp-ot's ferret_b9 ... ferret_b13
// (emp-ot docs/lpn-security.md), primal LPN with regular noise over a 10-local
// code, every counted attack at >= 141 bits over F_2 and over large fields.
// (t, logk, tree_depth): k = 2^logk, block size 2^tree_depth, N = t * 2^tree_depth.
//
//   set   table   t      logk  depth  N per round   lowest cost (bits)
//   b9    1 MiB   2,400  16    9       1,228,800    144.2
//   b10   2 MiB   2,800  17    10      2,867,200    141.3
//   b11   4 MiB   3,800  18    11      7,782,400    142.8
//   b12   8 MiB   4,200  19    12     17,203,200    141.2   (default)
//   b13  16 MiB   4,500  20    13     36,864,000    141.4
//
// Every set bootstraps through two ferret_b9 stages (pre0 from the COPE base
// sVOLE, then pre), which cover any set's t + k + 1 base pairs
// (<= 1,053,077 <= ferret_b9's N), as emp-ot's F_p sVOLE runs a ferret_b9
// pre-stage. Each stage and each round reads its own public matrix (lpn.h).
constexpr PrimalLPNParameterFp61 primal_lpn_set(std::size_t t, int logk, int depth) {
  return PrimalLPNParameterFp61(t << depth, t, std::size_t{1} << logk, depth,
                                2400u << 9, 2400, std::size_t{1} << 16, 9,
                                2400u << 9, 2400, std::size_t{1} << 16, 9);
}
const static PrimalLPNParameterFp61 ferret_b9  = primal_lpn_set(2400, 16, 9);
const static PrimalLPNParameterFp61 ferret_b10 = primal_lpn_set(2800, 17, 10);
const static PrimalLPNParameterFp61 ferret_b11 = primal_lpn_set(3800, 18, 11);
const static PrimalLPNParameterFp61 ferret_b12 = primal_lpn_set(4200, 19, 12);
const static PrimalLPNParameterFp61 ferret_b13 = primal_lpn_set(4500, 20, 13);
const static PrimalLPNParameterFp61 fp_default = ferret_b12;

// Domains of the per-stage public LPN matrices (lpn.h), shared with
// PrimalValueLocal (mvole.h), which regenerates a receiver's values.
inline constexpr char kLpnMain[] = "vole:lpn:main";
inline constexpr char kLpnPre[]  = "vole:lpn:pre";
inline constexpr char kLpnPre0[] = "vole:lpn:pre0";

template <typename IO, typename FP, typename FPS> 
class VoleTriple {
public:
  IO *io;
  IO **ios;
  int party;
  std::size_t threads;
  PrimalLPNParameterFp61 param;
  std::size_t noise_type;
  std::size_t M;
  std::size_t ot_used, ot_limit;
  bool is_malicious;
  FP *pre_yz_send = nullptr;
  FPS *pre_yz_recv = nullptr;

  ProgBaseCot<IO> *cot;
  OTPre<IO> *pre_ot = nullptr;

  FP Delta;
  LpnFp<10> *lpn = nullptr;
  ThreadPool *pool = nullptr;
  MpfssRegFp<IO, FP, FPS> *mpfss = nullptr;

  PRG prg;
  PRG prog_prg;

  // Public LPN-matrix root of this instance (lpn.h derives one matrix per
  // stage and round from it). The receiver (BOB) owns the LPN secret, so it
  // picks the root, uniformly unless set_lpn_seed() fixed it, and sends it to
  // ALICE at the start of setup(); choosing the matrix gives ALICE nothing.
  block lpn_root = zero_block;
  bool lpn_root_set = false;

  // Passive per-phase timing accumulators (microseconds) over all extend
  // rounds, for bench/: base COT for the GGM levels, MPFSS, LPN, tail copy.
  double t_cot = 0, t_mpfss = 0, t_lpn = 0, t_copy = 0, t_round = 0;

  VoleTriple(int party, std::size_t threads, IO **ios,
             PrimalLPNParameterFp61 param = fp_default) {
    this->io = ios[0];
    this->threads = threads;
    this->party = party;
    this->ios = ios;
    this->param = param;

    cot = new ProgBaseCot<IO>(party, io, true);
    cot->cot_gen_pre();

    pool = new ThreadPool(threads);
  }

  ~VoleTriple() {
    if (pre_yz_send != nullptr)
      delete[] pre_yz_send;
    if (pre_yz_recv != nullptr)
      delete[] pre_yz_recv;
    if (pre_ot != nullptr)
      delete pre_ot;
    if (lpn != nullptr)
      delete lpn;
    if (pool != nullptr)
      delete pool;
    if (mpfss != nullptr)
      delete mpfss;
    if (cot != nullptr)
      delete cot;
  }

  void setup_prog(block seed) {
    prog_prg.reseed(&seed);
  }

  // BOB only, before setup(): use this matrix root instead of a fresh one
  // (n-party VOLE, where one value vector must be reproducible elsewhere).
  void set_lpn_seed(block root) {
    if (party != BOB) error("set_lpn_seed: the receiver (BOB) picks the LPN matrix");
    lpn_root = root;
    lpn_root_set = true;
  }

  // Agree on the matrix root: BOB sends, ALICE receives.
  void exchange_lpn_root() {
    if (party == BOB) {
      if (!lpn_root_set) prg.random_block(&lpn_root, 1);
      io->send_data(&lpn_root, sizeof(block));
      io->flush();
    } else {
      io->recv_data(&lpn_root, sizeof(block));
    }
    lpn_root_set = true;
  }

  template<typename S>
  void setup(FP delta) {
    this->Delta = delta;
    setup<S>();
  }

  FP delta() {
    if (party == ALICE)
      return this->Delta;
    else {
      error("No delta for BOB");
      return 0;
    }
  }

  void extend_initialization() {
    lpn = new LpnFp<10>(param.n, param.k, pool, pool->size(), lpn_root, kLpnMain);
    mpfss = new MpfssRegFp<IO, FP, FPS>(party, threads, param.n, param.t,
                               param.log_bin_sz, pool, ios);
    mpfss->set_malicious();

    pre_ot = new OTPre<IO>(io, mpfss->tree_height - 1, mpfss->tree_n);
    // Base pairs reserved per round: the MPFSS point values (none for the
    // constant-1 bit instantiation) + the check mask + the LPN secret.
    M = mpfss->base_pairs() + param.k;
    ot_limit = param.n - M;
    ot_used = ot_limit;
  }

  // sender
  void extend_send(FP *buffer) {
    auto _t = clock_start();
    cot->prog_cot_gen(pre_ot, pre_ot->n);
    t_cot += time_from(_t); _t = clock_start();
    mpfss->sender_init(Delta);
    mpfss->mpfss_sender(buffer, pre_yz_send, pre_ot);
    t_mpfss += time_from(_t); _t = clock_start();
    lpn->compute(buffer, pre_yz_send + mpfss->base_pairs());
    t_lpn += time_from(_t); _t = clock_start();
    std::copy(buffer + ot_limit, buffer + ot_limit + M, pre_yz_send);
    t_copy += time_from(_t);
  }

  // receiver
  void extend_recv(FPS *buffer_mac) {
    auto _t = clock_start();
    bool *pre_bool_ini = new bool[pre_ot->n];
    prog_prg.random_bool(pre_bool_ini, pre_ot->n);
    cot->prog_cot_gen(pre_ot, pre_bool_ini, pre_ot->n);
    delete[] pre_bool_ini;
    t_cot += time_from(_t); _t = clock_start();

    mpfss->recver_init();
    mpfss->mpfss_recver(buffer_mac, pre_yz_recv, pre_ot);
    mpfss->set_vec_x(buffer_mac, pre_yz_recv);
    t_mpfss += time_from(_t); _t = clock_start();
    lpn->compute(buffer_mac, pre_yz_recv + mpfss->base_pairs());
    t_lpn += time_from(_t); _t = clock_start();
    std::copy(buffer_mac + ot_limit, buffer_mac + ot_limit + M, pre_yz_recv);
    t_copy += time_from(_t);
  }

  // store mac and val separately
  template<typename S>
  void setup() {
    exchange_lpn_root();
    extend_initialization();

    // pre-processing tools
    LpnFp<10> lpn_pre0(param.n_pre0, param.k_pre0, pool, pool->size(), lpn_root, kLpnPre0);
    MpfssRegFp<IO, FP, FPS> mpfss_pre0(party, threads, param.n_pre0, param.t_pre0,
                              param.log_bin_sz_pre0, pool, ios);
    mpfss_pre0.set_malicious();
    OTPre<IO> pre_ot_ini0(ios[0], mpfss_pre0.tree_height - 1,
                          mpfss_pre0.tree_n);

    // generate tree_n*(depth-1) COTs
    // generate 2*tree_n+k_pre triples and extend
    std::size_t M_pre0 = pre_ot_ini0.n;
    BaseSvoleFp<IO, FP> *svole0;
    std::size_t triple_n0 = mpfss_pre0.base_pairs() + param.k_pre0;
    FP *pre_yz0_send = nullptr;
    FPS *pre_yz0_recv = nullptr;
    if (party == ALICE) {
      cot->prog_cot_gen(&pre_ot_ini0, M_pre0);

      FP *key = new FP[triple_n0];
      svole0 = new BaseSvoleFp<IO, FP>(party, ios[0], Delta);
      svole0->compute_send64(key, triple_n0);
      vole_traits<FPS>::normalize_base_keys(key, triple_n0);

      pre_yz0_send = new FP[param.n_pre0];
      mpfss_pre0.sender_init(Delta);
      mpfss_pre0.mpfss_sender(pre_yz0_send, key, &pre_ot_ini0);
      lpn_pre0.compute(pre_yz0_send, key + mpfss_pre0.base_pairs());
      delete[] key;

    } else {

      bool *pre_bool_ini = new bool[pre_ot_ini0.n];
      prog_prg.random_bool(pre_bool_ini, pre_ot_ini0.n);
      cot->prog_cot_gen(&pre_ot_ini0, pre_bool_ini, M_pre0);
      delete[] pre_bool_ini;

      FPS *mac = new FPS[triple_n0];
      FP *x = new FP[triple_n0+1];
      svole0 = new BaseSvoleFp<IO, FP>(party, ios[0]);

      if(FP::PR_num_pack == 1) {
        for(std::size_t i = 0; i < triple_n0; ++i)
          vole_traits<FPS>::rand_base_value(x[i], prog_prg);
        // x[triple_n0] is the base-sVOLE check mask `a` (recver_check64): it
        // must be a uniform FIELD element even when the base values are bits,
        // or the sender, who picks chi, would learn sum(chi_i * x_i) up to
        // one bit, i.e. 128 linear relations on the receiver's base bits.
        x[triple_n0].rand(prog_prg);
      } else {
        S *buf = new S[triple_n0+1];
        for(std::size_t i = 0; i < triple_n0+1; ++i) {
          buf[i].rand(prog_prg);
          x[i].assign_no_mod(FP::copy_compose(buf[i].val));
        }
        delete[] buf;
      }

      svole0->template compute_recv64<FPS>(mac, x, triple_n0);

      pre_yz0_recv = new FPS[param.n_pre0];
      mpfss_pre0.recver_init();
      mpfss_pre0.mpfss_recver(pre_yz0_recv, mac, &pre_ot_ini0);
      mpfss_pre0.set_vec_x(pre_yz0_recv, mac);
      lpn_pre0.compute(pre_yz0_recv, mac + mpfss_pre0.base_pairs());
      delete[] mac;
      delete[] x;

    }
    delete svole0;

    // pre-processing tools
    LpnFp<10> lpn_pre(param.n_pre, param.k_pre, pool, pool->size(), lpn_root, kLpnPre);
    MpfssRegFp<IO, FP, FPS> mpfss_pre(party, threads, param.n_pre, param.t_pre,
                             param.log_bin_sz_pre, pool, ios);
    mpfss_pre.set_malicious();
    OTPre<IO> pre_ot_ini(ios[0], mpfss_pre.tree_height - 1, mpfss_pre.tree_n);

    // generate tree_n*(depth-1) COTs
    // generate 2*tree_n+k_pre triples and extend
    std::size_t M_pre = pre_ot_ini.n;
    if(party == ALICE) {
      cot->prog_cot_gen(&pre_ot_ini, M_pre);

      pre_yz_send = new FP[param.n_pre];
      mpfss_pre.sender_init(Delta);
      mpfss_pre.mpfss_sender(pre_yz_send, pre_yz0_send, &pre_ot_ini);
      lpn_pre.compute(pre_yz_send, pre_yz0_send + mpfss_pre.base_pairs());
      delete[] pre_yz0_send;
    } else {
      bool *pre_bool_ini = new bool[pre_ot_ini.n];
      prog_prg.random_bool(pre_bool_ini, pre_ot_ini.n);
      cot->prog_cot_gen(&pre_ot_ini, pre_bool_ini, M_pre);
      delete[] pre_bool_ini;

      pre_yz_recv = new FPS[param.n_pre];
      mpfss_pre.recver_init();
      mpfss_pre.mpfss_recver(pre_yz_recv, pre_yz0_recv, &pre_ot_ini);
      mpfss_pre.set_vec_x(pre_yz_recv, pre_yz0_recv);
      lpn_pre.compute(pre_yz_recv, pre_yz0_recv + mpfss_pre.base_pairs());
      delete[] pre_yz0_recv;
    }

  }

  std::size_t extend_inplace_send(FP *data_yz, std::size_t byte_space) {
    if (byte_space < param.n)
      error("space not enough");
    std::size_t tp_output_n = byte_space - M;
    if (tp_output_n % ot_limit != 0)
      error("call byte_memory_need_inplace \
				to get the correct length of memory space");
    std::size_t round = tp_output_n / ot_limit;
    FP *pt = data_yz;
    for (std::size_t i = 0; i < round; ++i) {
      auto _r = clock_start();
      extend_send(pt);
      t_round += time_from(_r);
      pt += ot_limit;
    }
    return tp_output_n;
  }

  std::size_t extend_inplace_recv(FPS *data_yz, std::size_t byte_space) {
    if (byte_space < param.n)
      error("space not enough");
    std::size_t tp_output_n = byte_space - M;
    if (tp_output_n % ot_limit != 0)
      error("call byte_memory_need_inplace \
				to get the correct length of memory space");
    std::size_t round = tp_output_n / ot_limit;
    FPS *pt = data_yz;
    for (std::size_t i = 0; i < round; ++i) {
      auto _r = clock_start();
      extend_recv(pt);
      t_round += time_from(_r);
      pt += ot_limit;
    }
    return tp_output_n;
  }

  std::size_t byte_memory_need_inplace(std::size_t tp_need) {
    std::size_t round = (tp_need - 1) / ot_limit;
    return round * ot_limit + param.n;
  }

  // debug function
  void check_triple(FP x, FP *y, std::size_t size, IO *io) {
    io->send_data(&x, sizeof(FP));
    io->send_data(y, size * sizeof(FP));
  }

  // debug function
  void check_triple(FPS *z, std::size_t size, IO *io) {
    FP delta;
    FP *k = new FP[size];
    io->recv_data(&delta, sizeof(FP));
    io->recv_data(k, size * sizeof(FP));
    for (std::size_t i = 0; i < size; ++i) {
      FP tmp = delta * z[i].getHigh();
      tmp = tmp + k[i];
      if (tmp != z[i].getLow()) {
        error("svole fails\n");
      }
    }

    delete[] k;
  }
};
#endif // _ITERATIVE_COT_H_
