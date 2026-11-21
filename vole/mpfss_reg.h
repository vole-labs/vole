#ifndef MPFSS_REG_FP_H__
#define MPFSS_REG_FP_H__

#include "vole/spfss_sender.h"
#include "vole/spfss_recver.h"
#include "vole/preot.h"
#include "emp-tool/emp-tool.h"
#include <set>
#include <algorithm>
#include <cstring>
#include <type_traits>

using namespace emp;

template <typename IO, typename FP, typename FPS>
class MpfssRegFp {
public:
  int party;
  std::size_t threads;
  std::size_t item_n, idx_max, m;
  std::size_t tree_height, leave_n;
  std::size_t tree_n;
  std::size_t nth;  // worker threads actually used = std::min(threads, tree_n) >= 1
  bool is_malicious;

  PRG prg;
  IO *netio;
  IO **ios;
  FP secret_share_x;
  FP *check_chialpha_buf = nullptr, *check_VW_buf = nullptr;
  FP *triple_yz_sender = nullptr;
  FPS *triple_yz_recver = nullptr;
  ThreadPool *pool;
  std::vector<uint32_t> item_pos_recver;

  MpfssRegFp(int party, std::size_t threads, std::size_t n, std::size_t t, std::size_t log_bin_sz,
             ThreadPool *pool, IO **ios) {
    this->party = party;
    this->threads = threads;
    this->netio = ios[0];
    this->ios = ios;

    this->pool = pool;
    this->is_malicious = false;

    // make sure n = t * leave_n
    this->item_n = t;
    this->idx_max = n;
    this->tree_height = log_bin_sz + 1;
    this->leave_n = 1 << (this->tree_height - 1);
    this->tree_n = this->item_n;
    // Never split tree_n trees over more than tree_n workers: with
    // threads > tree_n the old width = tree_n / threads was 0 and
    // ios[th] divided by zero.
    this->nth = std::max<std::size_t>(1, std::min(threads, tree_n));

    if (party == BOB)
      check_chialpha_buf = new FP[item_n];
    check_VW_buf = new FP[item_n];
  }

  ~MpfssRegFp() {
    if (check_chialpha_buf != nullptr)
      delete[] check_chialpha_buf;
    delete[] check_VW_buf;
  }

  void set_malicious() { is_malicious = true; }

  void sender_init(FP delta) { secret_share_x = delta; }

  void recver_init() { item_pos_recver.resize(this->item_n); }

  // Point values consumed from the base pairs (0 when the point value is the
  // constant 1, as in Ferret's SPCOT) and base pairs forming the check mask.
  static constexpr bool unit_point = vole_traits<FPS>::unit_point_value;
  static constexpr std::size_t point_pairs_per_tree = unit_point ? 0 : 1;
  static constexpr std::size_t mask_pairs = vole_traits<FPS>::mask_pairs;
  // Base pairs a caller must supply: point values (if any) + the mask.
  std::size_t base_pairs() const { return tree_n * point_pairs_per_tree + mask_pairs; }

  void set_vec_x(FPS *out, FPS *in) {
    uint32_t ptr = 0;
    for (std::size_t i = 0; i < tree_n; ++i) {
      uint32_t pt = ptr + item_pos_recver[i];
      out[pt].setHigh(unit_point ? FP(1) : in[i].getHigh());
      ptr += leave_n;
    }
  }

  // Per-tree "gamma" key / mask pair: the base pair's key/MAC for a field
  // value, or Delta / zero for the constant point value 1 (in characteristic
  // 2, secret_sum = Delta - sum(leaves) = Delta xor sum(leaves)).
  FP tree_gamma(const FP *triple_yz_key, std::size_t i) const {
    return unit_point ? secret_share_x : triple_yz_key[i];
  }
  FPS tree_delta2(const FPS *triple_yz_val_mac, std::size_t i) const {
    return unit_point ? FPS() : triple_yz_val_mac[i];
  }

  // sender
  void mpfss_sender(FP *sparse_vector,
             FP *triple_yz_key,
             OTPre<IO> *ot) {
    std::vector<SpfssSenderFp<IO, FP> *> senders;
    std::vector<std::future<void>> fut;
    // All tree seeds from one PRG call; one scratch tree per worker (only used
    // when FP is not block-sized, e.g. FP61's 8-byte keys).
    std::vector<block> seeds(tree_n);
    prg.random_block(seeds.data(), tree_n);
    std::size_t width = tree_n / nth;
    std::vector<block *> scratch(nth, nullptr);
    if (sizeof(FP) != sizeof(block))
      for (std::size_t th = 0; th < nth; ++th) scratch[th] = new block[leave_n];
    for (std::size_t i = 0; i < tree_n; ++i) {
      std::size_t th = (width == 0) ? 0 : std::min(i / width, nth - 1);
      senders.push_back(new SpfssSenderFp<IO, FP>(netio, tree_height, seeds[i], scratch[th]));
      ot->choices_sender();
    }
    netio->flush();
    ot->reset();

    std::size_t start = 0, end = width;
    for (std::size_t th = 0; th + 1 < nth; ++th) {
      fut.push_back(pool->enqueue(
          [this, th, start, end, senders, ot, sparse_vector, triple_yz_key]() {
            for (auto i = start; i < end; ++i) {
              senders[i]->compute(sparse_vector + i * leave_n, secret_share_x, tree_gamma(triple_yz_key, i));
              senders[i]->template send<OTPre<IO>>(ot, ios[th], i);
              ios[th]->flush();
            }
          }));
      start = end;
      end += width;
    }
    end = tree_n;
    for (auto i = start; i < end; ++i) {
      senders[i]->compute(sparse_vector + i * leave_n, secret_share_x, tree_gamma(triple_yz_key, i));
      senders[i]->template send<OTPre<IO>>(ot, ios[nth - 1], i);
      ios[nth - 1]->flush();
    }
    for (auto &f : fut)
      f.get();

    if (is_malicious) {
      // Ring (Z_{2^k}) only: GGM-tree tag check (MozZ2karella Fig. 5 step 6),
      // a no-op branch for the field types. Batched: 3 messages total for all
      // tree_n trees (K_top buffer -> shared challenge seed -> Gamma buffer)
      // instead of a per-tree ping-pong, so the round-trips are O(1), not
      // O(tree_n). The per-leaf tag/Gamma compute is still threaded.
      if (field_traits<FP>::is_ring) {
        block *ktop = new block[tree_n];
        ggm_parallel([&](std::size_t i) { ktop[i] = senders[i]->ggm_tag_top(); });
        netio->send_data(ktop, tree_n * sizeof(block));
        netio->flush();

        block chi_seed;
        netio->recv_data(&chi_seed, sizeof(block));
        block *xi = new block[leave_n];
        emp::uni_hash_coeff_gen(xi, chi_seed, (int)leave_n);

        block *gamma = new block[tree_n];
        ggm_parallel([&](std::size_t i) { gamma[i] = senders[i]->ggm_tag_response(xi); });
        netio->send_data(gamma, tree_n * sizeof(block));
        netio->flush();
        delete[] ktop;
        delete[] xi;
        delete[] gamma;
      }

      block *seed = new block[nth];
      seed_expand(seed, nth);
      std::vector<std::future<void>> fut;
      std::size_t start = 0, end = width;
      for (std::size_t th = 0; th + 1 < nth; ++th) {
        fut.push_back(
            pool->enqueue([this, th, start, end, senders, seed, sparse_vector]() {
              for (auto i = start; i < end; ++i) {
                check_VW_buf[i] = senders[i]->consistency_check_msg_gen(ios[th], seed[th], sparse_vector + i * leave_n, i);
              }
            }));
        start = end;
        end += width;
      }
      end = tree_n;
      for (auto i = start; i < end; ++i) {
        check_VW_buf[i] = senders[i]->consistency_check_msg_gen(ios[nth - 1], seed[nth - 1], sparse_vector + i * leave_n, i);
      }
      for (auto &f : fut)
        f.get();
      delete[] seed;
      consistency_batch_check(vole_traits<FPS>::pack_mask_key(triple_yz_key + tree_n * point_pairs_per_tree), tree_n);
    }

    for (auto p : senders)
      delete p;
    for (auto s : scratch)
      if (s != nullptr) delete[] s;
  }

  void mpfss_recver(FPS *sparse_vector,
             FPS *triple_yz_val_mac,
             OTPre<IO> *ot) {

    std::vector<SpfssRecverFp<IO, FP> *> recvers;
    std::vector<std::future<void>> fut;
    for (std::size_t i = 0; i < tree_n; ++i) {
      recvers.push_back(new SpfssRecverFp<IO, FP>(netio, tree_height));
      ot->choices_recver(recvers[i]->b);
      item_pos_recver[i] = recvers[i]->get_index();
    }
    netio->flush();
    ot->reset();

    std::size_t width = tree_n / nth;
    std::size_t start = 0, end = width;
    for (std::size_t th = 0; th + 1 < nth; ++th) {
      fut.push_back(pool->enqueue(
          [this, th, start, end, recvers, ot, sparse_vector, triple_yz_val_mac]() {
            for (auto i = start; i < end; ++i) {
              recvers[i]->template recv<OTPre<IO>>(ot, ios[th], i);
              recvers[i]->compute(sparse_vector + i * leave_n, tree_delta2(triple_yz_val_mac, i));
              ios[th]->flush();
            }
          }));
      start = end;
      end += width;
    }
    end = tree_n;
    for (auto i = start; i < end; ++i) {
      recvers[i]->template recv<OTPre<IO>>(ot, ios[nth - 1], i);
      recvers[i]->compute(sparse_vector + i * leave_n, tree_delta2(triple_yz_val_mac, i));
      ios[nth - 1]->flush();
    }
    for (auto &f : fut)
      f.get();

    if (is_malicious) {
      // Ring (Z_{2^k}) only: GGM-tree tag check (MozZ2karella Fig. 5 step 6),
      // a no-op branch for the field types. Batched into 3 messages total
      // (K_top buffer -> shared challenge seed -> Gamma buffer); abort if any
      // tree's Gamma mismatches.
      if (field_traits<FP>::is_ring) {
        block *ktop = new block[tree_n];
        netio->recv_data(ktop, tree_n * sizeof(block));
        ggm_parallel([&](std::size_t i) { recvers[i]->ggm_tag_top(ktop[i]); });

        block chi_seed;
        prg.random_block(&chi_seed, 1);
        netio->send_data(&chi_seed, sizeof(block));
        netio->flush();
        block *xi = new block[leave_n];
        emp::uni_hash_coeff_gen(xi, chi_seed, (int)leave_n);

        block *gamma_mine = new block[tree_n];
        ggm_parallel([&](std::size_t i) { gamma_mine[i] = recvers[i]->ggm_tag_response(xi); });
        block *gamma_recv = new block[tree_n];
        netio->recv_data(gamma_recv, tree_n * sizeof(block));

        bool ok = true;
        for (std::size_t i = 0; i < tree_n; ++i)
          if (!cmpBlock(&gamma_recv[i], &gamma_mine[i], 1)) ok = false;
        delete[] ktop;
        delete[] xi;
        delete[] gamma_mine;
        delete[] gamma_recv;
        if (!ok) error("MPFSS GGM-tree tag check fails");
      }

      block *seed = new block[nth];
      seed_expand(seed, nth);
      std::vector<std::future<void>> fut;
      std::size_t start = 0, end = width;
      for (std::size_t th = 0; th + 1 < nth; ++th) {
        fut.push_back(
            pool->enqueue([this, th, start, end, recvers, seed, sparse_vector, triple_yz_val_mac]() {
              for (auto i = start; i < end; ++i) {
                recvers[i]->consistency_check_msg_gen(
                    check_chialpha_buf[i], check_VW_buf[i],
                    ios[th], triple_yz_val_mac[i], seed[th], sparse_vector + i * leave_n, i);
              }
            }));
        start = end;
        end += width;
      }
      end = tree_n;
      for (auto i = start; i < end; ++i) {
        recvers[i]->consistency_check_msg_gen(
            check_chialpha_buf[i], check_VW_buf[i], ios[nth - 1],
            triple_yz_val_mac[i], seed[nth - 1], sparse_vector + i * leave_n, i);
      }
      for (auto &f : fut)
        f.get();
      delete[] seed;

      consistency_batch_check(triple_yz_val_mac,
          vole_traits<FPS>::pack_mask(triple_yz_val_mac + tree_n * point_pairs_per_tree), tree_n);
    }

    for (auto p : recvers)
      delete p;
  }

  // Run fn(i) for every tree i in [0, tree_n), split across the thread pool.
  // Used to parallelize the (I/O-free) GGM tag/Gamma compute before the batched
  // sends. Uses nth = std::min(threads, tree_n) workers.
  template <typename Fn>
  void ggm_parallel(Fn fn) {
    std::size_t width = tree_n / nth;
    std::vector<std::future<void>> fut;
    std::size_t s = 0, e = width;
    for (std::size_t th = 0; th + 1 < nth; ++th) {
      fut.push_back(pool->enqueue([s, e, fn]() {
        for (std::size_t i = s; i < e; ++i) fn(i);
      }));
      s = e;
      e += width;
    }
    for (std::size_t i = s; i < tree_n; ++i) fn(i);
    for (auto &f : fut) f.get();
  }

  void seed_expand(block *seed, std::size_t threads) {
    block sd = zero_block;
    if (party == ALICE) {
      netio->recv_data(&sd, sizeof(block));
    } else {
      prg.random_block(&sd, 1);
      netio->send_data(&sd, sizeof(block));
      netio->flush();
    }
    PRG prg2(&sd);
    prg2.random_block(seed, threads);
  }

  // Final check, Delta-holder side. vb = sum V_i - (y + Delta * x*) is linear
  // in Delta, so it is compared with the receiver's va by an equality test in
  // which this side commits first: commit to vb under a fresh 128-bit salt,
  // receive va, and open (send the salt) only if va == vb. A receiver that
  // deviated in x* holds only the hiding commitment, nothing to test
  // candidate Deltas against offline, and the opening never leaves on a
  // mismatch. (With a plain hash of vb, a receiver sending x* + e learned
  // H(va + e * Delta): a 2^61 offline search for FP61's Delta.) The
  // commitment binds vb before va arrives, so a sender that corrupted trees
  // cannot adapt to va; the va it then sees belongs to a round that aborts.
  void consistency_batch_check(FP y, std::size_t num) {
    FP x_star;
    netio->recv_data(&x_star, sizeof(FP));
    require_canonical(&x_star, 1, "MPFSS check: non-canonical value from the receiver");
    FP tmp = secret_share_x * x_star;
    tmp = y + tmp;
    FP vb = tmp.negate(); // y_star

    for (std::size_t i = 0; i < num; ++i)
      vb = vb + check_VW_buf[i];
    block salt, com[2];
    prg.random_block(&salt, 1);
    check_commit(vb, salt, com);
    if (test_cheat == kCheatBadOpening) salt = salt ^ makeBlock(0, 1);  // test hook
    netio->send_data(com, sizeof(com));
    netio->flush();
    FP va;
    netio->recv_data(&va, sizeof(FP));
    require_canonical(&va, 1, "MPFSS check: non-canonical value from the receiver");
    if (std::memcmp(&va, &vb, sizeof(FP)) != 0)
      error("MPFSS batch check fails");
    netio->send_data(&salt, sizeof(block));
    netio->flush();
  }

  // values are in delta2 (or the constant 1 for the bit instantiation)
  // z is the (packed) mask
  // Final check, receiver side (see the Delta-holder side above): send x*,
  // receive the commitment to vb, reveal va, and check the opening.
  template <typename MaskT>
  void consistency_batch_check(FPS *delta2, MaskT z, std::size_t num) {
    FP beta_mul_chialpha((uint64_t)0);
    for (std::size_t i = 0; i < num; ++i) {
      FP beta = unit_point ? FP(1) : delta2[i].getHigh();
      FP tmp = beta * check_chialpha_buf[i];
      beta_mul_chialpha = beta_mul_chialpha + tmp;
    }
    FP x_star = beta_mul_chialpha.negate();
    x_star = z.getHigh() + x_star;
    if (test_cheat == kCheatXStar) x_star = x_star + FP((uint64_t)1);  // test hook
    if constexpr (std::is_integral_v<decltype(x_star.val)>)               // test hook
      if (test_cheat == kCheatNonCanonical) x_star.val += (decltype(x_star.val))FP::PR;
    netio->send_data(&x_star, sizeof(FP));
    netio->flush();

    FP va = z.getLow();
    va = va.negate();
    for (std::size_t i = 0; i < num; ++i)
      va = va + check_VW_buf[i];
    if (test_cheat == kCheatVa) va = va + FP((uint64_t)1);  // test hook

    // The commitment arrives before va leaves; the sender opens it only if
    // the values match.
    block com[2];
    netio->recv_data(com, sizeof(com));
    netio->send_data(&va, sizeof(FP));
    netio->flush();
    block salt;
    netio->recv_data(&salt, sizeof(block));
    block opened[2];
    check_commit(va, salt, opened);
    if (!cmpBlock(com, opened, 2))
      error("MPFSS check: the sender's opening does not match its commitment");
  }

  // Commitment of the final check: RO(kDomCheck, salt, v), two blocks,
  // hiding while the 128-bit salt is secret.
  static void check_commit(const FP &v, block salt, block out[2]) {
    RO(kDomCheck, salt).absorb(&v, sizeof(FP)).squeeze_digest(out);
  }
  static constexpr char kDomCheck[] = "vole:mpfss:check";

  // Test hooks (test/mpfss_check.cpp); 0 = honest.
  int test_cheat = 0;
  static constexpr int kCheatXStar = 1;         // receiver: x* + 1
  static constexpr int kCheatNonCanonical = 2;  // receiver: x* + p (integer-stored fields)
  static constexpr int kCheatVa = 3;            // receiver: va + 1
  static constexpr int kCheatBadOpening = 4;    // sender: opens with another salt

};
#endif
