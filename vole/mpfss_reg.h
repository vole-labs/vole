#ifndef MPFSS_REG_FP_H__
#define MPFSS_REG_FP_H__

#include "vole/spfss_sender.h"
#include "vole/spfss_recver.h"
#include "vole/preot.h"
#include "emp-tool/emp-tool.h"
#include <set>
#include <algorithm>

using namespace emp;

template <typename IO, typename FP, typename FPS>
class MpfssRegFp {
public:
  int party;
  std::size_t threads;
  std::size_t item_n, idx_max, m;
  std::size_t tree_height, leave_n;
  std::size_t tree_n;
  std::size_t nth;  // worker threads actually used = min(threads, tree_n) >= 1
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

  void set_vec_x(FPS *out, FPS *in) {
    uint32_t ptr = 0;
    for (std::size_t i = 0; i < tree_n; ++i) {
      uint32_t pt = ptr + item_pos_recver[i];
      out[pt].setHigh(in[i].getHigh());
      ptr += leave_n;
    }
  }

  // sender
  void mpfss_sender(FP *sparse_vector,
             FP *triple_yz_key,
             OTPre<IO> *ot) {
    vector<SpfssSenderFp<IO, FP> *> senders;
    vector<future<void>> fut;
    for (std::size_t i = 0; i < tree_n; ++i) {
      senders.push_back(new SpfssSenderFp<IO, FP>(netio, tree_height));
      ot->choices_sender();
    }
    netio->flush();
    ot->reset();

    std::size_t width = tree_n / nth;
    std::size_t start = 0, end = width;
    for (std::size_t th = 0; th + 1 < nth; ++th) {
      fut.push_back(pool->enqueue(
          [this, th, start, end, senders, ot, sparse_vector, triple_yz_key]() {
            for (auto i = start; i < end; ++i) {
              senders[i]->compute(sparse_vector + i * leave_n, secret_share_x, triple_yz_key[i]);
              senders[i]->template send<OTPre<IO>>(ot, ios[th], i);
              ios[th]->flush();
            }
          }));
      start = end;
      end += width;
    }
    end = tree_n;
    for (auto i = start; i < end; ++i) {
      senders[i]->compute(sparse_vector + i * leave_n, secret_share_x, triple_yz_key[i]);
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
      vector<future<void>> fut;
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
      consistency_batch_check(triple_yz_key[tree_n], tree_n);
    }

    for (auto p : senders)
      delete p;
  }

  void mpfss_recver(FPS *sparse_vector,
             FPS *triple_yz_val_mac,
             OTPre<IO> *ot) {

    vector<SpfssRecverFp<IO, FP> *> recvers;
    vector<future<void>> fut;
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
              recvers[i]->compute(sparse_vector + i * leave_n, triple_yz_val_mac[i]);
              ios[th]->flush();
            }
          }));
      start = end;
      end += width;
    }
    end = tree_n;
    for (auto i = start; i < end; ++i) {
      recvers[i]->template recv<OTPre<IO>>(ot, ios[nth - 1], i);
      recvers[i]->compute(sparse_vector + i * leave_n, triple_yz_val_mac[i]);
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
      vector<future<void>> fut;
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

      consistency_batch_check(triple_yz_val_mac, triple_yz_val_mac[tree_n], tree_n);
    }

    for (auto p : recvers)
      delete p;
  }

  // Run fn(i) for every tree i in [0, tree_n), split across the thread pool.
  // Used to parallelize the (I/O-free) GGM tag/Gamma compute before the batched
  // sends. Uses nth = min(threads, tree_n) workers.
  template <typename Fn>
  void ggm_parallel(Fn fn) {
    std::size_t width = tree_n / nth;
    vector<future<void>> fut;
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

  void consistency_batch_check(FP y, std::size_t num) {
    FP x_star;
    x_star.recv(netio);
    FP tmp = secret_share_x * x_star;
    tmp = y + tmp;
    FP vb = tmp.negate(); // y_star

    for (std::size_t i = 0; i < num; ++i)
      vb = vb + check_VW_buf[i];
    block h = vb.hash();
    netio->send_data(&h, sizeof(block));
    netio->flush();
  }

  // values are in delta2
  // z is the mask
  void consistency_batch_check(FPS *delta2, FPS z, std::size_t num) {
    FP beta_mul_chialpha((uint64_t)0);
    for (std::size_t i = 0; i < num; ++i) {
      FP tmp = delta2[i].getHigh() * check_chialpha_buf[i];
      beta_mul_chialpha = beta_mul_chialpha + tmp;
    }
    FP x_star = beta_mul_chialpha.negate();
    x_star = z.getHigh() + x_star;
    x_star.template send<IO>(netio);
    netio->flush();

    FP va = z.getLow();
    va = va.negate();
    for (std::size_t i = 0; i < num; ++i)
      va = va + check_VW_buf[i];

    block h = va.hash();
    block r;
    netio->recv_data(&r, sizeof(block));
    if (!cmpBlock(&r, &h, 1))
      error("MPFSS batch check fails");
  }

};
#endif
