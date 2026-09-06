#ifndef PROG_COT_H__
#define PROG_COT_H__

#include "emp-ot/emp-ot.h"

template<typename IO>
class ProgBaseCot : public BaseCot<IO> { public:
  using BaseCot<IO>::party;
	using BaseCot<IO>::one;
  using BaseCot<IO>::minusone;
	using BaseCot<IO>::ot_delta;
	using BaseCot<IO>::io;
	using BaseCot<IO>::iknp;
	using BaseCot<IO>::malicious;

	ProgBaseCot(int party, IO *io, bool malicious = false) :
      BaseCot<IO>(party, io, malicious){

	}
	
	~ProgBaseCot() {
	
  }

	void prog_cot_gen(OTPre<IO> *pre_ot, std::size_t size) {
		block *ot_data = new block[size];
    iknp->send_cot(ot_data, size);
    io->flush();
    for(std::size_t i = 0; i < size; ++i)
      ot_data[i] = ot_data[i] & minusone;
    pre_ot->send_pre(ot_data, ot_delta);
		delete[] ot_data;
	}

	void prog_cot_gen(OTPre<IO> *pre_ot, bool *pre_bool_ini, std::size_t size) {
		block *ot_data = new block[size];
    PRG prg;
    iknp->recv_cot(ot_data, pre_bool_ini, size);
    block ch[2];
    ch[0] = zero_block;
    ch[1] = makeBlock(0, 1);
    for(std::size_t i = 0; i < size; ++i)
      ot_data[i] = 
          (ot_data[i] & minusone) ^ ch[pre_bool_ini[i]];
    pre_ot->recv_pre(ot_data, pre_bool_ini);
		delete[] ot_data;
	}


	// debug
	bool check_cot(block *data, std::size_t len) {
		if(party == ALICE) {
			io->send_block(&ot_delta, 1);
			io->send_block(data, len); 
			io->flush();
			return true;
		} else {
			block * tmp = new block[len];
			block ch[2];
			io->recv_block(ch+1, 1);
			ch[0] = zero_block;
			io->recv_block(tmp, len);
			for(std::size_t i = 0; i < len; ++i)
				tmp[i] = tmp[i] ^ ch[getLSB(data[i])];
			bool res = cmpBlock(tmp, data, len);
			delete[] tmp;
			return res;
		}
	}
};

#endif
