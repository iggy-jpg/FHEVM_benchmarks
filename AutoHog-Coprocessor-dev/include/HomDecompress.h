#pragma once

#include <tfhe++.hpp>
#include <vector>
#include <chrono>
#include <omp.h>
#include "drift_ms_tfhepp.hpp" 
#include "drift_ms_params.hpp" 
#include "Autohog_Cop.hpp"  
#include "HomPack.h"        // TLWE2TRLWEPacking
using namespace TFHEpp;

namespace HomDecompress {

using Clock       = std::chrono::high_resolution_clock;
using MilliDouble = std::chrono::duration<double, std::milli>;

using DriftKeyLvl0 = driftms::DriftMsKey<lvl0param>;
using DriftKeyLvl2 = driftms::DriftMsKey<lvl2param>;


using Pin  = TFHEpp::lvl1param;
using Pout = TFHEpp::lvl2param;

// ===================================================================

// ===================================================================
inline double Decompress_PBSiks_lvl1_to_lvl2(
    TFHEpp::TRLWE<Pout>                                &packed_out,   
    const TFHEpp::TRLWE<Pin>                           &packed_in,    
    const TFHEpp::KeySwitchingKey<TFHEpp::lvl10param>  &ksk10,        
    const TFHEpp::BootstrappingKeyFFT<TFHEpp::lvl02param> &bkfft,     // 0 -> 2 PBS
    const TFHEpp::AnnihilateKey<Pout>                  &ahk2,         
    int                                                 slot_count)   
{
    using namespace TFHEpp;

    
    std::vector<TLWE<Pin>> slots_l1 = TRLWE_Decompress<Pin>(packed_in, slot_count);
    const int S = slot_count;

    std::vector<TLWE<lvl0param>> slots_l0(S);
    std::vector<TLWE<Pout>>      slots_l2(S);

    auto t0 = Clock::now();
    
    Polynomial<Pout> poly;
    for (auto &co : poly) co = Pout::μ;
    
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < S; ++i) {
        // 2.1 IKS: lvl1 -> lvl0
        IdentityKeySwitch<lvl10param>(slots_l0[i], slots_l1[i], ksk10);

        // 2.2 PBS: lvl0 -> lvl2

    
        TRLWE<lvl02param::targetP> acc;
        BlindRotate<lvl02param>(acc, slots_l0[i], bkfft, poly);

        SampleExtractIndex<lvl02param::targetP>(slots_l2[i], acc, 0);
    }

    
    HomPack::TLWE2TRLWEPacking<Pout>(packed_out, slots_l2, ahk2);

    auto t1 = Clock::now();
    return MilliDouble(t1 - t0).count();
}

// ===================================================================


// ===================================================================
inline double Decompress_PBS_lvl1_to_lvl2(
    TFHEpp::TRLWE<Pout>                                &packed_out,   
    const TFHEpp::TRLWE<Pin>                           &packed_in,    
    const TFHEpp::BootstrappingKeyFFT<TFHEpp::lvl12param> &bkfft12,   // 1 -> 2 PBS
    const TFHEpp::AnnihilateKey<Pout>                  &ahk2,         
    int                                                 slot_count)   
{
    using namespace TFHEpp;

    
    std::vector<TLWE<Pin>> slots_l1 = TRLWE_Decompress<Pin>(packed_in, slot_count);
    const int S = slot_count;

    std::vector<TLWE<Pout>> slots_l2(S);

    auto t0 = Clock::now();
    
    Polynomial<Pout> poly;
    for (auto &co : poly) co = Pout::μ;
    
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < S; ++i) {
        
        

        
        TRLWE<lvl12param::targetP> acc;
        BlindRotate<lvl12param>(acc, slots_l1[i], bkfft12, poly);

        
        SampleExtractIndex<lvl12param::targetP>(slots_l2[i], acc, 0);
    }

    
    HomPack::TLWE2TRLWEPacking<Pout>(packed_out, slots_l2, ahk2);

    auto t1 = Clock::now();
    return MilliDouble(t1 - t0).count();
}


inline double Decompress_DirectIKS_lvl1_to_lvl2(
    TFHEpp::TRLWE<Pout>                              &packed_out,
    const TFHEpp::TRLWE<Pin>                         &packed_in,
    const TFHEpp::KeySwitchingKey<TFHEpp::lvl12param> &ksk12,
    const TFHEpp::AnnihilateKey<Pout>                &ahk2,
    int                                                slot_count,
    const DriftKeyLvl2 *drift_key_lvl2 = nullptr)
{
    using namespace TFHEpp;

    std::vector<TLWE<Pin>> slots_l1 = TRLWE_Decompress<Pin>(packed_in, slot_count);
    const int S = slot_count;

    std::vector<TLWE<Pout>> slots_l2(S);

    auto t0 = Clock::now();

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < S; ++i) {
        IdentityKeySwitch<lvl12param>(slots_l2[i], slots_l1[i], ksk12);
        //centered::apply_centered_correction_to_body<Pout>(slots_l2[i], /*log_mod=*/12);
        // if (drift_key_lvl2 != nullptr) {
        
        //     //printf("Drifting slot %d at lvl2...\n", i);
        //     // driftms::improve_lwe_ciphertext_modulus_switch_noise_for_binary_key<Pout>(
        //     //     slots_l2[i], *drift_key_lvl2
        //     // );
        //     centered::apply_centered_correction_to_body<Pout>(slots_l2[i], drift_key_lvl2->log_mod);
        // }
    }

    HomPack::TLWE2TRLWEPacking<Pout>(packed_out, slots_l2, ahk2);

    auto t1 = Clock::now();
    return MilliDouble(t1 - t0).count();
}

} // namespace HomDecompress
