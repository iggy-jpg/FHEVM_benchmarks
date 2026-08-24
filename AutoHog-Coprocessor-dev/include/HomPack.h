#pragma once
#include <vector>
#include <tfhe++.hpp>

namespace HomPack {


template <class P>
inline void TLWE2TRLWEPacking(TFHEpp::TRLWE<P>& out,
                              std::vector<TFHEpp::TLWE<P>>& in,
                              const TFHEpp::AnnihilateKey<P>& ahk)
{
    TFHEpp::TLWE2TRLWEPacking<P>(out, in, ahk);
}


template <class P>
inline void TLWE2TRLWEPacking(TFHEpp::TRLWE<P>& out,
                              const std::vector<TFHEpp::TLWE<P>>& in_const,
                              const TFHEpp::AnnihilateKey<P>& ahk)
{
    auto in = in_const; 
    TFHEpp::TLWE2TRLWEPacking<P>(out, in, ahk);
}

} // namespace HomPack

template <class P>
inline std::vector<TFHEpp::TLWE<P>> TRLWE_Decompress(const TFHEpp::TRLWE<P>& ct, int count) {
    std::vector<TFHEpp::TLWE<P>> result(count);
    for (int i = 0; i < count; i++) {
        TFHEpp::SampleExtractIndex<P>(result[i], ct, i);
    }
    return result;
}