#ifndef DIST_PSI_CVOLE_FP_H__
#define DIST_PSI_CVOLE_FP_H__

// Field-generic *committed* VOLE.
//
// Dual-LPN (expand-accumulate) VOLE templated on a finite field or ring (see
// vole/fields/), plus the LPN commitment of the paper (Sec. 4.1). Built from
// the field-generic silent-VOLE components MpfssRegFp (regular multi-point
// FSS), BaseSvoleFp (COPE-based base sVOLE), ProgBaseCot/OTPre, and the
// expand-accumulate code AccumFp + LpnFpEA.
//
// Why dual-LPN: a committed VOLE must publish a commitment that *binds* the
// VOLE value x. With dual LPN, x = H * e_u is fully determined by the regular
// sparse seed e_u, so committing to e_u binds x. (A commitment over a dense x
// is a linear map with trivial collisions -> not binding; binding holds only
// over the low-weight e_u.) Unlike a primal VOLE there is no "A*seed" term, so
// the MPFSS only needs M = (t+1)+(t_com+1) base sVOLEs per extend.
//
// Base sVOLE bootstrap (Ferret-style self-seeding, like VoleTriple): setup()
// generates the first M base sVOLEs once via COPE; each extend seeds its two
// MPFSS from those, and reserves the tail M of its own VOLE output back into the
// seed buffer for the next extend. So COPE runs exactly once (in setup), not
// per extend.
//
// Roles (matching MpfssRegFp / BaseSvoleFp):
//   ALICE = verifier   : holds Delta and the keys K[x], K[com].
//   BOB   = committer   : holds the values x, com and the macs M[x], M[com],
//                         such that M = K + value * Delta.
//
// Commitment (paper Sec. 4.1, dual-LPN commitment): a single shared buffer
// holds the RAW regular-sparse [e_u || e_r]. The main output is
// x = H * accumulate(e_u) (expand-accumulate dual LPN), and the commitment is
//   com = H_u * e_u + H_r * e_r
// computed by ComMatrixFp directly on the raw sparse vector (no accumulator on
// the commitment path): H_u is column-sparse random (every column hit), H_r is
// dense random (the paper's H_2). The same authenticated e_u feeds both, so com
// is consistent with the VOLE on x for free. The consistency check verifies
// Hash(M[com]) == Hash(K[com] + com * Delta).
//
// Binding (paper Thm. 1): a difference of two openings has at most 2(t + t_com)
// nonzeros and 2^{2 t log_bin + 2 t_com log_bin_com} possible supports, so
//   n_com  >=  2(t + t_com) + (2 t log_bin + 2 t_com log_bin_com + lambda) / log|F|
// is required (lambda = statistical parameter). Hiding: H_r * e_r is a dual-LPN
// sample of a random dense n_com x N_com matrix, the best linear test succeeds
// with probability about (n_com / N_com)^{t_com}.

#include "vole/vole-fp.h"
#include "vole/base_cot.h"
#include "vole/preot.h"
#include "vole/accumulator_fp.h"
#include "vole/lpn_ea.h"
#include "vole/com_matrix.h"
#include <algorithm>
#include <cmath>

using namespace emp;

// Fixed, distinct PRP seeds for the three public matrices (H must not share a
// seed with the commitment matrices, or their rows become correlated).
const static block cvole_seed_lpn   = makeBlock(0x564f4c452d4c504eULL, 1);  // "VOLE-LPN"
const static block cvole_seed_com_u = makeBlock(0x434f4d2d48555f55ULL, 1);  // "COM-HU"
const static block cvole_seed_com_r = makeBlock(0x434f4d2d48525f52ULL, 1);  // "COM-HR"

class CVoleFpParam {
public:
  // main VOLE: output length n (=|x|), t nonzeros, leave size 2^log_bin.
  // e_u length N = t * 2^log_bin (>= n, dual/compressing expansion).
  std::size_t n, t, log_bin;
  // commitment: output length n_com (=|com|), t_com nonzeros, leave 2^log_bin_com.
  // e_r length N_com = t_com * 2^log_bin_com. n_com == 0 means "derive from the
  // field and the binding bound" (see cvole_n_com_for / cvole_resolve_param).
  std::size_t n_com, t_com, log_bin_com;
  // commitment matrix H_u: entries per column (ignored if dense_com). A
  // structural dependency needs col_weight + 1 columns whose supports all lie
  // inside a common col_weight-row set (probability 1 / C(n_com, col_weight)
  // per column, ~2^-150 at 32 rows of 832), and the generic threshold is
  // log2(n_com) ~ 10; 32 keeps a 3x margin and costs one PRP block plus one
  // multiply-add per entry. Both parties must use the same value.
  std::size_t col_weight;
  bool dense_com;

  CVoleFpParam() {}
  CVoleFpParam(std::size_t n, std::size_t t, std::size_t log_bin,
               std::size_t n_com, std::size_t t_com, std::size_t log_bin_com,
               std::size_t col_weight = 32, bool dense_com = false)
      : n(n), t(t), log_bin(log_bin),
        n_com(n_com), t_com(t_com), log_bin_com(log_bin_com),
        col_weight(col_weight), dense_com(dense_com) {}

  std::size_t N() const { return t << log_bin; }            // |e_u|
  std::size_t N_com() const { return t_com << log_bin_com; } // |e_r|
};

// ---- LPN parameters ----
// Main expansion (n, t, log_bin): output length n = 2^XX, regular weight
// t = 224, log_bin = XX - 5, so |e_u| = 224 * 2^(XX-5) = 7 * n (code rate 1/7).
// These are the paper's conservative EA-code parameters (Table 1: w = 224,
// L = 15 / 19 / 23 for 2^20 / 2^24 / 2^28), which match the C = 10, C*ln(N)
// taps-per-row expander implemented by LpnFpEA. (The previous t = 160 at rate
// 1/5 gave only ~75 bits by the paper's Appendix A.1 estimate.)
//
// Commitment (t_com, log_bin_com) = (128, 4), |e_r| = 2048, n_com = 0 = AUTO:
//   n_com is derived per FIELD and per scale at construction time from the
//   binding bound above (cvole_n_com_for<FP>), rounded up to a multiple of 64:
//     2^20 : 768 (F_2^128), 832 (FP107), 832 (FP61)
//     2^28 : 832 (F_2^128), 832 (FP107), 896 (FP61)
//   hiding  : (n_com / 2048)^128 <= 2^-152 for the H_r * e_r mask.
// (The old (160, 100, 3) commitment had n_com < t + t_com and was not binding.)
const static CVoleFpParam cvole_fp_n2to20 = CVoleFpParam(1ull << 20, 224, 15, 0, 128, 4);
const static CVoleFpParam cvole_fp_n2to21 = CVoleFpParam(1ull << 21, 224, 16, 0, 128, 4);
const static CVoleFpParam cvole_fp_n2to22 = CVoleFpParam(1ull << 22, 224, 17, 0, 128, 4);
const static CVoleFpParam cvole_fp_n2to23 = CVoleFpParam(1ull << 23, 224, 18, 0, 128, 4);
const static CVoleFpParam cvole_fp_n2to24 = CVoleFpParam(1ull << 24, 224, 19, 0, 128, 4);
const static CVoleFpParam cvole_fp_n2to25 = CVoleFpParam(1ull << 25, 224, 20, 0, 128, 4);
const static CVoleFpParam cvole_fp_n2to26 = CVoleFpParam(1ull << 26, 224, 21, 0, 128, 4);
const static CVoleFpParam cvole_fp_n2to27 = CVoleFpParam(1ull << 27, 224, 22, 0, 128, 4);
const static CVoleFpParam cvole_fp_n2to28 = CVoleFpParam(1ull << 28, 224, 23, 0, 128, 4);

// Setup-scale params: small + fast, for tests (same rate 1/7 and weight 224, so
// N = 7 * 2^14). The commitment parameters are the real ones.
const static CVoleFpParam cvole_fp_n2to14 = CVoleFpParam(1ull << 14, 224, 9, 0, 128, 4);

const static CVoleFpParam cvole_fp_default = cvole_fp_n2to20;

// (#5) Pick the smallest tabulated param whose per-round size n covers n_need
// (mirrors find_vole_scale in psi.h / nvole.h). Beyond 2^28 it returns the
// largest entry and the driver falls back to multiple rounds.
inline CVoleFpParam cvole_fp_param_for(std::size_t n_need) {
  if (n_need <= (1ull << 14)) return cvole_fp_n2to14;  // small/test scale
  if (n_need <= (1ull << 20)) return cvole_fp_n2to20;
  if (n_need <= (1ull << 21)) return cvole_fp_n2to21;
  if (n_need <= (1ull << 22)) return cvole_fp_n2to22;
  if (n_need <= (1ull << 23)) return cvole_fp_n2to23;
  if (n_need <= (1ull << 24)) return cvole_fp_n2to24;
  if (n_need <= (1ull << 25)) return cvole_fp_n2to25;
  if (n_need <= (1ull << 26)) return cvole_fp_n2to26;
  if (n_need <= (1ull << 27)) return cvole_fp_n2to27;
  return cvole_fp_n2to28;
}

// Smallest n_com satisfying the Theorem-1 binding bound for field FP (log|F| =
// FP::PR_bit_len), statistical parameter lambda, rounded up to a multiple of 64.
// Over the ring Z_{2^k} this is only a placeholder (the bound is a field
// argument).
template <typename FP>
inline std::size_t cvole_n_com_for(const CVoleFpParam &p, std::size_t lambda = 40) {
  double bits = 2.0 * p.t * p.log_bin + 2.0 * p.t_com * p.log_bin_com + (double)lambda;
  std::size_t need = 2 * (p.t + p.t_com) +
                     (std::size_t)std::ceil(bits / (double)FP::PR_bit_len);
  return ((need + 63) / 64) * 64;
}

// Fill in n_com when the table entry leaves it at 0 (auto). Every party must
// resolve with the same FP, which the drivers guarantee by sharing the field.
template <typename FP>
inline CVoleFpParam cvole_resolve_param(CVoleFpParam p) {
  if (p.n_com == 0) p.n_com = cvole_n_com_for<FP>(p);
  return p;
}

template <typename IO, typename FP, typename FPS>
class CVoleFp {
public:
  int party;
  std::size_t threads;
  IO **ios;
  IO *io;
  CVoleFpParam param;

  FP Delta;
  PRG prog_prg;

  ThreadPool *pool = nullptr;
  ProgBaseCot<IO> *cot = nullptr;
  BaseSvoleFp<IO, FP> *svole = nullptr;

  MpfssRegFp<IO, FP, FPS> *mpfss = nullptr;      // e_u
  MpfssRegFp<IO, FP, FPS> *mpfss_com = nullptr;  // e_r
  OTPre<IO> *ot_pre = nullptr;
  OTPre<IO> *ot_pre_com = nullptr;

  // Expand-accumulate (EA) code for x (accumulate e_u, then sparse expand), and
  // the commitment matrix applied to the raw [e_u || e_r].
  AccumFp *acc_com = nullptr;         // accumulate e_u in place (length N)
  LpnFpEA *lpn = nullptr;             // H   : A(e_u) (N)         -> x   (n)
  ComMatrixFp<FP> *hcom = nullptr;    // [H_u | H_r] : [e_u||e_r] -> com (n_com)

  // Debug/test hook: the sparse positions (into [e_u || e_r]) and values used by
  // the most recent extend_local round.
  std::vector<std::size_t> dbg_pos;
  std::vector<FP> dbg_val;

  // Reserved base sVOLE seed, bootstrapped once in setup() and refreshed from
  // each extend's output tail (Ferret-style self-seeding).
  FP *pre_yz_send = nullptr;   // ALICE: M base-sVOLE keys
  FPS *pre_yz_recv = nullptr;  // BOB:   M base-sVOLE (value||mac)
  std::size_t cu = 0, cr = 0;  // base sVOLEs for the e_u / e_r MPFSS (tree_n+1)
  std::size_t M = 0;           // = cu + cr, reserved & refreshed each extend
  std::size_t ot_limit = 0;    // usable VOLE outputs per extend = n - M

  // Passive per-phase timing accumulators (microseconds), summed across all
  // extend rounds. Used by the bench/ harness to split runtime into the VOLE
  // (value) pipeline, the commitment pipeline, and the shared accumulate. The
  // overhead is a handful of clock reads per round (negligible vs the work), so
  // they are always on and harmless to the correctness tests in test/.
  double t_vole = 0;   // A: cot + MPFSS(e_u) + LPN H  -> x
  double t_com = 0;    // B: cot + MPFSS(e_r) + [H_u|H_r] -> com
  double t_accum = 0;  // accumulate over e_u

  // King (local) mode: reproduces the committer's VALUE side (u and com) from
  // the seed alone, with no IO / COPE / COT. The
  // king chose all seeds, so it can publish each party's commitment locally.
  bool is_local = false;
  FP *pre_yz_val_local = nullptr;  // king: M base-sVOLE values, self-seeded

  CVoleFp(int party, std::size_t threads, IO **ios,
          CVoleFpParam param = cvole_fp_default)
      : party(party), threads(threads), ios(ios),
        param(cvole_resolve_param<FP>(param)) {
    this->io = ios[0];
    pool = new ThreadPool(threads);
    cot = new ProgBaseCot<IO>(party, io, true);
    cot->cot_gen_pre();
  }

  // Local (non-interactive) constructor for the king. No IO/COPE/COT: it only
  // reproduces the committer's value side from prog_prg via setup_local /
  // extend_local.
  CVoleFp(std::size_t threads, CVoleFpParam param = cvole_fp_default)
      : party(-1), threads(threads), ios(nullptr), io(nullptr),
        param(cvole_resolve_param<FP>(param)), is_local(true) {
    pool = new ThreadPool(threads);
  }

  ~CVoleFp() {
    if (svole != nullptr) delete svole;
    if (mpfss != nullptr) delete mpfss;
    if (mpfss_com != nullptr) delete mpfss_com;
    if (ot_pre != nullptr) delete ot_pre;
    if (ot_pre_com != nullptr) delete ot_pre_com;
    if (acc_com != nullptr) delete acc_com;
    if (lpn != nullptr) delete lpn;
    if (hcom != nullptr) delete hcom;
    if (pre_yz_send != nullptr) delete[] pre_yz_send;
    if (pre_yz_recv != nullptr) delete[] pre_yz_recv;
    if (pre_yz_val_local != nullptr) delete[] pre_yz_val_local;
    if (cot != nullptr) delete cot;
    if (pool != nullptr) delete pool;
  }

  void setup_prog(block seed) { prog_prg.reseed(&seed); }

  // Draw n field values from prog_prg (the value side of the base sVOLE). Used
  // by both the committer's recv_base_svole_prog and the king's setup_local so
  // they consume prog_prg identically.
  void gen_prog_values(FP *xin, std::size_t n) {
    if (FP::PR_num_pack == 1) {
      for (std::size_t i = 0; i < n; ++i)
        xin[i].rand(prog_prg);
    } else {
      uint64_t *buf = new uint64_t[n];
      prog_prg.random_data_unaligned(buf, n * sizeof(uint64_t));
      for (std::size_t i = 0; i < n; ++i) {
        buf[i] = buf[i] & FP::PR_mask;
        xin[i].assign_no_mod(FP::copy_compose(buf[i]));
      }
      delete[] buf;
    }
  }

  // Committer-side base sVOLE whose VALUES are derived from prog_prg, so that
  // the same prog_seed reproduces the same e_u (hence the same x and com)
  // across sessions / verifiers. Mirrors BaseSvoleFp's no-x overload but feeds
  // prog_prg instead of a fresh PRG. `cnt` correlations need cnt+1 inputs.
  void recv_base_svole_prog(FPS *mac, std::size_t cnt) {
    FP *xin = new FP[cnt + 1];
    gen_prog_values(xin, cnt + 1);  // values; macs come from COPE below
    svole->template compute_recv64<FPS>(mac, xin, cnt);
    delete[] xin;
  }

  void setup(FP delta) {
    this->Delta = delta;
    setup();
  }

  // The three public matrices (shared by the interactive and king-local paths).
  void init_matrices() {
    acc_com = new AccumFp(param.N(), pool, pool->size());
    lpn = new LpnFpEA(param.n, param.N(), pool, pool->size(), cvole_seed_lpn);
    hcom = new ComMatrixFp<FP>(param.n_com, param.N(), param.N_com(),
                               param.col_weight, cvole_seed_com_u,
                               cvole_seed_com_r, pool, pool->size(),
                               param.dense_com);
  }

  void setup() {
    mpfss = new MpfssRegFp<IO, FP, FPS>(party, threads, param.N(), param.t,
                                        param.log_bin, pool, ios);
    mpfss->set_malicious();
    mpfss_com = new MpfssRegFp<IO, FP, FPS>(party, threads, param.N_com(),
                                            param.t_com, param.log_bin_com,
                                            pool, ios);
    mpfss_com->set_malicious();

    ot_pre = new OTPre<IO>(io, mpfss->tree_height - 1, mpfss->tree_n);
    ot_pre_com = new OTPre<IO>(io, mpfss_com->tree_height - 1, mpfss_com->tree_n);

    init_matrices();

    if (party == ALICE)
      svole = new BaseSvoleFp<IO, FP>(party, io, Delta);
    else
      svole = new BaseSvoleFp<IO, FP>(party, io);

    cu = mpfss->tree_n + 1;       // base sVOLEs for the e_u MPFSS
    cr = mpfss_com->tree_n + 1;   // base sVOLEs for the e_r MPFSS
    M = cu + cr;
    ot_limit = param.n - M;

    // One-time base sVOLE bootstrap (the only COPE call). Subsequent extends
    // self-seed from the reserved tail of their own output. For BOB the values
    // come from prog_prg, so the committed input is reproducible across fresh
    // instances seeded with the same prog_seed (call setup_prog before setup).
    if (party == ALICE) {
      pre_yz_send = new FP[M];
      svole->compute_send64(pre_yz_send, M);
    } else {
      pre_yz_recv = new FPS[M];
      recv_base_svole_prog(pre_yz_recv, M);
    }
  }

  // ALICE (verifier): produce one round of key K[x] (length n) and
  // K[com] (length n_com).
  void extend_send(FP *x_key, FP *com_key) {
    std::size_t N = param.N(), N_com = param.N_com();
    FP *sparse = new FP[N + N_com];  // [e_u || e_r]  (vole_f2k's pre_sparse_yz)

    // MPFSS seeds come from the reserved base sVOLE (pre_yz_send), sliced as in
    // vole_f2k (pre_yz for e_u, pre_yz+t+1 for e_r) -- no per-extend COPE.
    auto _t = clock_start();
    cot->prog_cot_gen(ot_pre, ot_pre->n);
    mpfss->sender_init(Delta);
    mpfss->mpfss_sender(sparse, pre_yz_send, ot_pre);
    t_vole += time_from(_t);

    _t = clock_start();
    cot->prog_cot_gen(ot_pre_com, ot_pre_com->n);
    mpfss_com->sender_init(Delta);
    mpfss_com->mpfss_sender(sparse + N, pre_yz_send + cu, ot_pre_com);
    t_com += time_from(_t);

    // ----- commitment on the RAW sparse vector, before the accumulate -----
    _t = clock_start();
    hcom->compute(com_key, sparse);
    t_com += time_from(_t);

    // ----- EA expand for x: accumulate e_u in place, then sparse expand -----
    _t = clock_start();
    acc_com->compute(sparse);
    t_accum += time_from(_t);

    _t = clock_start();
    for (std::size_t i = 0; i < param.n; ++i) x_key[i].setZero();
    lpn->compute(x_key, sparse);
    t_vole += time_from(_t);

    // Reserve the tail M outputs as the next extend's base sVOLE seed.
    std::copy(x_key + ot_limit, x_key + ot_limit + M, pre_yz_send);

    delete[] sparse;
  }

  // BOB (committer): produce one round of x (val||mac, length n) and
  // com (val||mac, n_com). The value must always be materialized (set_vec_x),
  // since the reserved tail's value seeds the next extend.
  void extend_recv(FPS *x_vm, FPS *com_vm) {
    std::size_t N = param.N(), N_com = param.N_com();
    FPS *sparse = new FPS[N + N_com];  // [e_u || e_r]

    // MPFSS seeds come from the reserved base sVOLE (pre_yz_recv), sliced for
    // e_u / e_r. OT choices (positions) still come from prog_prg.
    auto _t = clock_start();
    {
      bool *pb = new bool[ot_pre->n];
      prog_prg.random_bool(pb, ot_pre->n);
      cot->prog_cot_gen(ot_pre, pb, ot_pre->n);
      delete[] pb;
    }
    mpfss->recver_init();
    mpfss->mpfss_recver(sparse, pre_yz_recv, ot_pre);
    mpfss->set_vec_x(sparse, pre_yz_recv);
    t_vole += time_from(_t);

    _t = clock_start();
    {
      bool *pb = new bool[ot_pre_com->n];
      prog_prg.random_bool(pb, ot_pre_com->n);
      cot->prog_cot_gen(ot_pre_com, pb, ot_pre_com->n);
      delete[] pb;
    }
    mpfss_com->recver_init();
    mpfss_com->mpfss_recver(sparse + N, pre_yz_recv + cu, ot_pre_com);
    mpfss_com->set_vec_x(sparse + N, pre_yz_recv + cu);
    t_com += time_from(_t);

    // ----- commitment on the RAW sparse vector, before the accumulate -----
    _t = clock_start();
    hcom->compute(com_vm, sparse);
    t_com += time_from(_t);

    // ----- EA expand for x (see extend_send) -----
    _t = clock_start();
    acc_com->compute(sparse);
    t_accum += time_from(_t);

    _t = clock_start();
    for (std::size_t i = 0; i < param.n; ++i) x_vm[i].setZero();
    lpn->compute(x_vm, sparse);
    t_vole += time_from(_t);

    // Reserve the tail M outputs as the next extend's base sVOLE seed.
    std::copy(x_vm + ot_limit, x_vm + ot_limit + M, pre_yz_recv);

    delete[] sparse;
  }

  // ---- king-side local computation (no IO/COPE/COT) ----
  // Reproduces the committer's VALUE side (u and com) from the seed. The king
  // chose the seeds, so it can compute every party's
  // commitment locally and publish it.
  void setup_local() {
    cu = param.t + 1;
    cr = param.t_com + 1;
    M = cu + cr;
    ot_limit = param.n - M;
    init_matrices();

    // Bootstrap the value chain, matching the committer's recv_base_svole_prog(M)
    // (which draws M+1 prog_prg values; the first M are the base sVOLE values).
    pre_yz_val_local = new FP[M];
    FP *xin = new FP[M + 1];
    gen_prog_values(xin, M + 1);
    for (std::size_t i = 0; i < M; ++i) pre_yz_val_local[i] = xin[i];
    delete[] xin;
  }

  // Place `values[j]` at the regular-sparse position of block j, where the
  // position is decoded from the OT choice bits exactly as SpfssRecverFp does.
  // `base` is the offset of this region inside [e_u || e_r]; the absolute
  // positions/values are appended to dbg_pos/dbg_val.
  void place_sparse(FP *out, const bool *pb, std::size_t num_trees,
                    std::size_t depth_minus1, std::size_t leave_n,
                    const FP *values, std::size_t base) {
    std::size_t ptr = 0;
    for (std::size_t j = 0; j < num_trees; ++j) {
      std::size_t pos = 0;
      for (std::size_t l = 0; l < depth_minus1; ++l) {
        pos <<= 1;
        if (!pb[ptr + l]) pos += 1;
      }
      out[j * leave_n + pos] = values[j];
      dbg_pos.push_back(base + j * leave_n + pos);
      dbg_val.push_back(values[j]);
      ptr += depth_minus1;
    }
  }

  // King: reproduce one round's value side x = u (length n) and com (n_com),
  // consuming prog_prg exactly as the committer's extend_recv (OT bits for e_u
  // then e_r) and self-seeding the same way.
  void extend_local(FP *x_val, FP *com_val) {
    std::size_t N = param.N(), N_com = param.N_com();
    std::size_t leave_u = (std::size_t)1 << param.log_bin;
    std::size_t leave_r = (std::size_t)1 << param.log_bin_com;
    FP *sparse = new FP[N + N_com];
    for (std::size_t i = 0; i < N + N_com; ++i) sparse[i].setZero();
    dbg_pos.clear();
    dbg_val.clear();

    // e_u value: positions from prog_prg OT bits, values = pre_yz_val_local[0:t]
    {
      std::size_t nb = param.log_bin * param.t;  // == ot_pre->n
      bool *pb = new bool[nb];
      prog_prg.random_bool(pb, nb);
      place_sparse(sparse, pb, param.t, param.log_bin, leave_u, pre_yz_val_local, 0);
      delete[] pb;
    }
    // e_r value: values = pre_yz_val_local[cu : cu+t_com]
    {
      std::size_t nb = param.log_bin_com * param.t_com;  // == ot_pre_com->n
      bool *pb = new bool[nb];
      prog_prg.random_bool(pb, nb);
      place_sparse(sparse + N, pb, param.t_com, param.log_bin_com, leave_r,
                   pre_yz_val_local + cu, N);
      delete[] pb;
    }

    // Commitment straight from the sparse (position, value) list.
    hcom->compute_sparse(com_val, dbg_pos.data(), dbg_val.data(), dbg_pos.size());

    acc_com->compute(sparse);
    for (std::size_t i = 0; i < param.n; ++i) x_val[i].setZero();
    lpn->compute(x_val, sparse);

    std::copy(x_val + ot_limit, x_val + ot_limit + M, pre_yz_val_local);
    delete[] sparse;
  }

  std::size_t extend_inplace_local(FP *x_val, FP *com_val, std::size_t rounds) {
    FP *pt = x_val;
    FP *pt_com = com_val;
    for (std::size_t r = 0; r < rounds; ++r) {
      extend_local(pt, pt_com);
      pt += ot_limit;
      pt_com += param.n_com;
    }
    scrub_tail(x_val, rounds);
    return x_usable(rounds);
  }

  // ---- iterated extension (Ferret-style self-seeding, like VoleTriple) ----
  // Each round yields ot_limit = n - M usable VOLE outputs (the tail M are
  // reserved as the next round's seed) plus n_com commitment outputs. For
  // `rounds` rounds the x buffer overlaps by M between consecutive rounds, so it
  // holds x_buf_size = (rounds-1)*ot_limit + n entries, of which only the first
  // x_usable = rounds*ot_limit are VOLE outputs for the caller. The last M
  // entries are the seed of the NEXT extend: extend_inplace_* copies them into
  // the internal seed buffer and then zeroes them in the caller's buffer, so
  // they can neither be consumed as correlations nor leak the next round's
  // seed. The com blocks are non-overlapping.
  std::size_t rounds_for(std::size_t n_need) const {
    return (n_need + ot_limit - 1) / ot_limit;
  }
  std::size_t x_buf_size(std::size_t rounds) const {
    return (rounds - 1) * ot_limit + param.n;
  }
  std::size_t x_usable(std::size_t rounds) const {
    return rounds * ot_limit;
  }
  std::size_t com_buf_size(std::size_t rounds) const {
    return rounds * param.n_com;
  }

  template <typename T>
  void scrub_tail(T *x, std::size_t rounds) {
    for (std::size_t i = x_usable(rounds); i < x_buf_size(rounds); ++i)
      x[i].setZero();
  }

  // Return the number of usable VOLE outputs (x_usable(rounds)).
  std::size_t extend_inplace_send(FP *x_key, FP *com_key, std::size_t rounds) {
    FP *pt = x_key;
    FP *pt_com = com_key;
    for (std::size_t r = 0; r < rounds; ++r) {
      extend_send(pt, pt_com);
      pt += ot_limit;
      pt_com += param.n_com;
    }
    scrub_tail(x_key, rounds);
    return x_usable(rounds);
  }

  std::size_t extend_inplace_recv(FPS *x_vm, FPS *com_vm, std::size_t rounds) {
    FPS *pt = x_vm;
    FPS *pt_com = com_vm;
    for (std::size_t r = 0; r < rounds; ++r) {
      extend_recv(pt, pt_com);
      pt += ot_limit;
      pt_com += param.n_com;
    }
    scrub_tail(x_vm, rounds);
    return x_usable(rounds);
  }

  // Reduce received field elements to their canonical encoding. A malicious
  // committer could otherwise send e.g. a >= p representative over FP61, which
  // is the same element but feeds mult_mod a 64-bit operand and would let two
  // clients hold byte-different copies of the same commitment. FP(val) always
  // reduces (a no-op for F_2^128 / Z_2^k).
  void canonicalize(FP *v, std::size_t n) {
    for (std::size_t i = 0; i < n; ++i) v[i] = FP(v[i].val);
  }

  // ---- committed-VOLE consistency check ----
  // BOB publishes the com VALUE (com_val, length com_len) and Hash(M[com]);
  // ALICE verifies Hash(K[com] + com*Delta) == Hash(M[com]). com_len = rounds *
  // n_com. com_val is passed in (rather than read from com_vm) so the committer
  // can publish one commitment value and reuse it with every client.
  void consistency_check_commit(const FP *com_val, FPS *com_vm,
                                std::size_t com_len) {  // BOB
    FP *com_mac = new FP[com_len];
    for (std::size_t i = 0; i < com_len; ++i)
      com_mac[i] = com_vm[i].getLow();
    block h[2];  // SHA-256 digest is 32 bytes
    Hash::hash_once(h, com_mac, com_len * sizeof(FP));
    io->send_data(com_val, com_len * sizeof(FP));
    io->send_data(h, 2 * sizeof(block));
    io->flush();
    delete[] com_mac;
  }

  void consistency_check_verify(FP *com_key, std::size_t com_len) {  // ALICE
    FP *com_val = new FP[com_len];
    block h_recv[2];
    io->recv_data(com_val, com_len * sizeof(FP));
    io->recv_data(h_recv, 2 * sizeof(block));
    canonicalize(com_val, com_len);

    FP *chk = new FP[com_len];
    for (std::size_t i = 0; i < com_len; ++i)
      chk[i] = com_key[i] + com_val[i] * Delta;
    block h[2];
    Hash::hash_once(h, chk, com_len * sizeof(FP));
    delete[] com_val;
    delete[] chk;
    if (memcmp(h, h_recv, 2 * sizeof(block)) != 0)
      error("CVoleFp consistency check failed");
  }

  // ---- consistency check against a PUBLISHED commitment (n-party / king) ----
  // Here the com VALUE is published by the king, not sent by the committer.
  // Committer (BOB) sends only Hash(M[com]); verifier (ALICE) checks
  // Hash(K[com] + com_pub*Delta) == Hash(M[com]) using the king's com_pub.
  void consistency_send_machash(FPS *com_vm, std::size_t com_len) {  // BOB
    FP *com_mac = new FP[com_len];
    for (std::size_t i = 0; i < com_len; ++i)
      com_mac[i] = com_vm[i].getLow();
    block h[2];
    Hash::hash_once(h, com_mac, com_len * sizeof(FP));
    io->send_data(h, 2 * sizeof(block));
    io->flush();
    delete[] com_mac;
  }

  void consistency_check_pub(FP *com_key, const FP *com_pub,
                             std::size_t com_len) {  // ALICE
    block h_recv[2];
    io->recv_data(h_recv, 2 * sizeof(block));
    FP *chk = new FP[com_len];
    for (std::size_t i = 0; i < com_len; ++i)
      chk[i] = com_key[i] + FP(com_pub[i].val) * Delta;  // canonical encoding
    block h[2];
    Hash::hash_once(h, chk, com_len * sizeof(FP));
    delete[] chk;
    if (memcmp(h, h_recv, 2 * sizeof(block)) != 0)
      error("CVoleFp published-commitment check failed");
  }

  // ---- debug: directly verify a VOLE relation M = K + val*Delta ----
  void check_vole_send(FP *key, std::size_t num) {  // ALICE
    io->send_data(&Delta, sizeof(FP));
    io->send_data(key, num * sizeof(FP));
    io->flush();
  }

  void check_vole_recv(FPS *vm, std::size_t num) {  // BOB
    FP delta;
    FP *key = new FP[num];
    io->recv_data(&delta, sizeof(FP));
    io->recv_data(key, num * sizeof(FP));
    for (std::size_t i = 0; i < num; ++i) {
      FP tmp = vm[i].getHigh() * delta;
      tmp = tmp + key[i];
      if (tmp != vm[i].getLow())
        error("CVoleFp VOLE relation check failed");
    }
    delete[] key;
  }
};

#endif // DIST_PSI_CVOLE_FP_H__
