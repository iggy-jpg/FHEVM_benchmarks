#pragma once
#include "HomGate.h"
#include <cstdlib>
#include <optional>
#include <string>
#include <vector>
#include "HomPKE.h"

namespace TFHEpp {

extern TLWE<lvl2param> CT_ZERO;
extern TLWE<lvl2param> CT_ONE;

template <typename TLWE_TYPE>
    inline void sanitize_tlwe(TLWE_TYPE& ct) {
        // Keep the historical workaround available, but allow consensus tests
        // to prove that the FFT path is deterministic without discarding bits.
        static const bool disabled =
            std::getenv("AUTOHOG_DISABLE_LOW4_SANITIZE") != nullptr;
        if (disabled) return;
        for (auto& val : ct) {
            val &= ~uint64_t{0xF};
        }
    }


void init_constants(const SecretKey& sk);
void init_constants(const PKE::SeededRlwePkPair<lvl2param>& spk);

struct HomAddResult {
    std::vector<TLWE<lvl2param>> SUM;
    std::optional<TLWE<lvl2param>> ovf_ct;
};

struct HomDivResult {
    std::vector<TLWE<lvl2param>> quotient;
    std::vector<TLWE<lvl2param>> remainder;
};

HomAddResult HomAdd(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<TLWE<lvl2param>>& B_in,
    int W,
    const EvalKey& ek,
    bool compute_overflow
);

HomAddResult HomSub(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<TLWE<lvl2param>>& B_in,
    int W,
    const EvalKey& ek,
    bool compute_overflow
);

//========================================================

TLWE<lvl2param> HomGT(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<TLWE<lvl2param>>& B_in,
    int W,
    const EvalKey& ek
);

TLWE<lvl2param> HomLE(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<TLWE<lvl2param>>& B_in,
    int W,
    const EvalKey& ek
);

TLWE<lvl2param> HomLT(
    const std::vector<TLWE<lvl2param>>& A_in,     
    const std::vector<TLWE<lvl2param>>& B_in,       
    int W,
    const EvalKey& ek
);

TLWE<lvl2param> HomGE(
    const std::vector<TLWE<lvl2param>>& A_in,     
    const std::vector<TLWE<lvl2param>>& B_in,       
    int W,
    const EvalKey& ek
);

TLWE<lvl2param> HomEQ(
    const std::vector<TLWE<lvl2param>>& A_in,     
    const std::vector<TLWE<lvl2param>>& B_in,       
    int W,
    const EvalKey& ek
);

TLWE<lvl2param> HomNE(
    const std::vector<TLWE<lvl2param>>& A_in,     
    const std::vector<TLWE<lvl2param>>& B_in,       
    int W,
    const EvalKey& ek
);

//========================================================

std::vector<TLWE<lvl2param>> HomAND(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<TLWE<lvl2param>>& B_in,
    int W,
    const EvalKey& ek,
    int threads_hint = -1);

std::vector<TLWE<lvl2param>> HomOR(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<TLWE<lvl2param>>& B_in,
    int W,
    const EvalKey& ek,
    int threads_hint = -1);

std::vector<TLWE<lvl2param>> HomXOR(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<TLWE<lvl2param>>& B_in,
    int W,
    const EvalKey& ek,
    int threads_hint = -1);

std::vector<TLWE<lvl2param>> HomNOT(
    const std::vector<TLWE<lvl2param>>& A_in,
    int W,
    const EvalKey& ek);

std::vector<TLWE<lvl2param>> HomSelect(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<TLWE<lvl2param>>& B_in,
    const TLWE<lvl2param>& sel, 
    int W,
    const EvalKey& ek,
    int threads_hint = -1);

std::vector<TLWE<lvl2param>> HomMul(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<TLWE<lvl2param>>& B_in,
    int W,
    const EvalKey& ek);

std::vector<TLWE<lvl2param>> HomMul_pal(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<TLWE<lvl2param>>& B_in,
    int W,
    const EvalKey& ek,
    int threads_hint = -1);

HomDivResult HomDiv(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<TLWE<lvl2param>>& B_in,
    int W,
    const EvalKey& ek);

HomDivResult HomDiv_pal(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<TLWE<lvl2param>>& B_in,
    int W,
    const EvalKey& ek,
    int threads_hint = -1);

 




HomDivResult HomDivNewton(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<TLWE<lvl2param>>& B_in,
    int W,
    const EvalKey& ek,
    int threads_hint = -1);

 



std::vector<TLWE<lvl2param>> HomColumnsCompressAdd(
    std::vector<std::vector<TLWE<lvl2param>>>& columns,
    int outW,
    const EvalKey& ek,
    bool parallel,
    int threads_hint = -1);

 




std::vector<TLWE<lvl2param>> HomMulPlain(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<uint8_t>& k_bits,
    int W,
    const EvalKey& ek);

std::vector<TLWE<lvl2param>> HomMulPlain_pal(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<uint8_t>& k_bits,
    int W,
    const EvalKey& ek,
    int threads_hint = -1);

 





HomDivResult HomDivPlain(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<uint8_t>& d_bits,
    int W,
    const EvalKey& ek);

HomDivResult HomDivPlain_pal(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<uint8_t>& d_bits,
    int W,
    const EvalKey& ek,
    int threads_hint = -1);

 




std::vector<std::vector<TLWE<lvl2param>>> HomMulPlainBatch(
    const std::vector<std::vector<TLWE<lvl2param>>>& As,
    const std::vector<uint8_t>& k_bits,
    int W,
    const EvalKey& ek,
    int threads_hint = -1);

 




std::vector<TLWE<lvl2param>> HomDivPlainQuotient(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<uint8_t>& d_bits,
    int W,
    const EvalKey& ek);

std::vector<TLWE<lvl2param>> HomDivPlainQuotient_pal(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<uint8_t>& d_bits,
    int W,
    const EvalKey& ek,
    int threads_hint = -1);

 


std::vector<std::vector<TLWE<lvl2param>>> HomDivPlainQuotientBatch(
    const std::vector<std::vector<TLWE<lvl2param>>>& As,
    const std::vector<uint8_t>& d_bits,
    int W,
    const EvalKey& ek,
    int threads_hint = -1);

 


std::vector<HomDivResult> HomDivPlainBatch(
    const std::vector<std::vector<TLWE<lvl2param>>>& As,
    const std::vector<uint8_t>& d_bits,
    int W,
    const EvalKey& ek,
    int threads_hint = -1);

 



TLWE<lvl2param> HomEQPlain(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<uint8_t>& k_bits,
    int W,
    const EvalKey& ek);

TLWE<lvl2param> HomEQPlain_pal(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<uint8_t>& k_bits,
    int W,
    const EvalKey& ek,
    int threads_hint = -1);

 


TLWE<lvl2param> HomGTPlain(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<uint8_t>& k_bits,
    int W,
    const EvalKey& ek);

TLWE<lvl2param> HomGTPlain_pal(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<uint8_t>& k_bits,
    int W,
    const EvalKey& ek,
    int threads_hint = -1);
TLWE<lvl2param> HomGEPlain(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<uint8_t>& k_bits,
    int W,
    const EvalKey& ek);

TLWE<lvl2param> HomGEPlain_pal(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<uint8_t>& k_bits,
    int W,
    const EvalKey& ek,
    int threads_hint = -1);
TLWE<lvl2param> HomLTPlain(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<uint8_t>& k_bits,
    int W,
    const EvalKey& ek);

TLWE<lvl2param> HomLTPlain_pal(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<uint8_t>& k_bits,
    int W,
    const EvalKey& ek,
    int threads_hint = -1);
TLWE<lvl2param> HomLEPlain(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<uint8_t>& k_bits,
    int W,
    const EvalKey& ek);

TLWE<lvl2param> HomLEPlain_pal(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<uint8_t>& k_bits,
    int W,
    const EvalKey& ek,
    int threads_hint = -1);

 


std::vector<TLWE<lvl2param>> HomANDPlain(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<uint8_t>& k_bits,
    int W);
std::vector<TLWE<lvl2param>> HomORPlain(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<uint8_t>& k_bits,
    int W);
std::vector<TLWE<lvl2param>> HomXORPlain(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<uint8_t>& k_bits,
    int W);

 



std::vector<TLWE<lvl2param>> HomSelectPlain(
    const std::vector<uint8_t>& a_bits,
    const std::vector<uint8_t>& b_bits,
    const TLWE<lvl2param>& sel,
    int W);


 



std::vector<TLWE<lvl2param>> HomMulPlainLow(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<uint8_t>& k_bits,
    int W,
    const EvalKey& ek);

std::vector<TLWE<lvl2param>> HomMulPlainLow_pal(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<uint8_t>& k_bits,
    int W,
    const EvalKey& ek,
    int threads_hint = -1);

 


std::vector<std::vector<TLWE<lvl2param>>> HomMulPlainLowBatch(
    const std::vector<std::vector<TLWE<lvl2param>>>& As,
    const std::vector<uint8_t>& k_bits,
    int W,
    const EvalKey& ek,
    int threads_hint = -1);

 


std::vector<TLWE<lvl2param>> HomMulCompressorLow(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<TLWE<lvl2param>>& B_in,
    int W,
    const EvalKey& ek);

std::vector<TLWE<lvl2param>> HomMulCompressorLow_pal(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<TLWE<lvl2param>>& B_in,
    int W,
    const EvalKey& ek,
    int threads_hint = -1);

 




std::vector<TLWE<lvl2param>> HomMulRectWindow(
    const std::vector<TLWE<lvl2param>>& A_in,
    int Wa,
    const std::vector<TLWE<lvl2param>>& B_in,
    int Wb,
    int lo_col,
    int hi_col,
    const EvalKey& ek,
    int threads_hint = -1);

std::vector<TLWE<lvl2param>> HomCast(
    const std::vector<TLWE<lvl2param>>& A_in,
    int fromW,
    int toW,
    bool signed_extend = false);

//========================================================

 




HomAddResult HomAdd_unsigned_wide(const std::vector<TLWE<lvl2param>>& A_in,
                                  const std::vector<TLWE<lvl2param>>& B_in,
                                  int W, const EvalKey& ek, bool compute_carry = false,
                                  int threads_hint = -1);
HomAddResult HomSub_unsigned_wide(const std::vector<TLWE<lvl2param>>& A_in,
                                  const std::vector<TLWE<lvl2param>>& B_in,
                                  int W, const EvalKey& ek, bool compute_borrow = false,
                                  int threads_hint = -1);
HomAddResult HomAdd_wide(const std::vector<TLWE<lvl2param>>& A_in,
                         const std::vector<TLWE<lvl2param>>& B_in,
                         int W, const EvalKey& ek, bool compute_overflow,
                         int threads_hint = -1);
HomAddResult HomSub_wide(const std::vector<TLWE<lvl2param>>& A_in,
                         const std::vector<TLWE<lvl2param>>& B_in,
                         int W, const EvalKey& ek, bool compute_overflow,
                         int threads_hint = -1);
TLWE<lvl2param> HomGT_unsigned_wide(const std::vector<TLWE<lvl2param>>& A_in,
                                    const std::vector<TLWE<lvl2param>>& B_in,
                                    int W, const EvalKey& ek, int threads_hint = -1);

HomAddResult HomAdd_pal(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<TLWE<lvl2param>>& B_in,
    int W,
    const EvalKey& ek,
    bool compute_overflow,
    int threads_hint = -1
);

HomAddResult HomSub_pal(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<TLWE<lvl2param>>& B_in,
    int W,
    const EvalKey& ek,
    bool compute_overflow,
    int threads_hint = -1
);

TLWE<lvl2param> HomGT_pal(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<TLWE<lvl2param>>& B_in,
    int W,
    const EvalKey& ek,
    int threads_hint = -1);

TLWE<lvl2param> HomLE_pal(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<TLWE<lvl2param>>& B_in,
    int W,
    const EvalKey& ek,
    int threads_hint = -1);

TLWE<lvl2param> HomLT_pal(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<TLWE<lvl2param>>& B_in,
    int W,
    const EvalKey& ek,
    int threads_hint = -1);

TLWE<lvl2param> HomGE_pal(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<TLWE<lvl2param>>& B_in,
    int W,
    const EvalKey& ek,
    int threads_hint = -1
);

TLWE<lvl2param> HomEQ_pal(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<TLWE<lvl2param>>& B_in,
    int W,
    const EvalKey& ek,
    int threads_hint = -1
);

TLWE<lvl2param> HomNE_pal(
    const std::vector<TLWE<lvl2param>>& A_in,
    const std::vector<TLWE<lvl2param>>& B_in,
    int W,
    const EvalKey& ek,
    int threads_hint = -1
);
 



HomAddResult HomAdd_unsigned(const std::vector<TLWE<lvl2param>> &A_in,
                                const std::vector<TLWE<lvl2param>> &B_in,
                                int W,
                                const EvalKey &ek,
                                bool compute_carry = false);

 



HomAddResult HomSub_unsigned(const std::vector<TLWE<lvl2param>> &A_in,
                                const std::vector<TLWE<lvl2param>> &B_in,
                                int W,
                                const EvalKey &ek,
                                bool compute_borrow = false);

 


HomAddResult HomAdd_unsigned_pal(const std::vector<TLWE<lvl2param>> &A_in,
                                    const std::vector<TLWE<lvl2param>> &B_in,
                                    int W,
                                    const EvalKey &ek,
                                    bool compute_carry = false,
                                    int threads_hint = -1);

 


HomAddResult HomSub_unsigned_pal(const std::vector<TLWE<lvl2param>> &A_in,
                                    const std::vector<TLWE<lvl2param>> &B_in,
                                    int W,
                                    const EvalKey &ek,
                                    bool compute_borrow = false,
                                    int threads_hint = -1);


TLWE<lvl2param> HomGT_unsigned(const std::vector<TLWE<lvl2param>> &A, const std::vector<TLWE<lvl2param>> &B, int W, const EvalKey &ek);

TLWE<lvl2param> HomLT_unsigned(const std::vector<TLWE<lvl2param>> &A, const std::vector<TLWE<lvl2param>> &B, int W, const EvalKey &ek);

TLWE<lvl2param> HomGE_unsigned(const std::vector<TLWE<lvl2param>> &A, const std::vector<TLWE<lvl2param>> &B, int W, const EvalKey &ek);

TLWE<lvl2param> HomLE_unsigned(const std::vector<TLWE<lvl2param>> &A, const std::vector<TLWE<lvl2param>> &B, int W, const EvalKey &ek);


TLWE<lvl2param> HomGT_unsigned_pal(const std::vector<TLWE<lvl2param>> &A, const std::vector<TLWE<lvl2param>> &B, int W, const EvalKey &ek, int threads_hint = -1);
TLWE<lvl2param> HomLT_unsigned_pal(const std::vector<TLWE<lvl2param>> &A, const std::vector<TLWE<lvl2param>> &B, int W, const EvalKey &ek, int threads_hint = -1);
TLWE<lvl2param> HomGE_unsigned_pal(const std::vector<TLWE<lvl2param>> &A, const std::vector<TLWE<lvl2param>> &B, int W, const EvalKey &ek, int threads_hint = -1);
TLWE<lvl2param> HomLE_unsigned_pal(const std::vector<TLWE<lvl2param>> &A, const std::vector<TLWE<lvl2param>> &B, int W, const EvalKey &ek, int threads_hint = -1);

}
