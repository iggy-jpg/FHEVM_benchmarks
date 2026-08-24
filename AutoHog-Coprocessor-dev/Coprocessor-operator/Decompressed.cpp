#include <tfhe++.hpp>
#include <vector>
#include <random>
#include <chrono>
#include <cstdio>
#include <omp.h>
#include <cstdlib>

#include "Autohog_Cop.hpp"
#include "HomPack.h"        // TRLWE_Decompress / TLWE2TRLWEPacking / annihilatekeygen
#include "HomPKE.h"         
#include "HomDecompress.h"  // Decompress_PBSiks_lvl1_to_lvl2 / Decompress_PBS_lvl1_to_lvl2 / Decompress_DirectIKS_lvl1_to_lvl2
#include "SecureRandom.h"

using namespace TFHEpp;

using Clock       = std::chrono::high_resolution_clock;
using MilliDouble = std::chrono::duration<double, std::milli>;


using Pin  = lvl1param;
using Pout = lvl2param;

int main(int argc, char **argv)
{
    
    omp_set_dynamic(0);
    omp_set_num_threads(omp_get_max_threads());

    const int NUM_TESTS = 100; 

    
    const std::vector<int> slot_list = {1, 32, 64, 128, 256, 512};

    // =====================
    
    // =====================
    SecretKey sk;

    EvalKey ek;
    ek.emplacebkfft<lvl02param>(sk); // PBS: lvl0 -> lvl2
    ek.emplacebkfft<lvl12param>(sk); 
    ek.emplaceiksk<lvl10param>(sk);  // IKS: lvl1 -> lvl0
    ek.emplaceiksk<lvl12param>(sk);  

    auto ksk10_sp = ek.get<KeySwitchingKey<lvl10param>>();
    auto &ksk10 = *ksk10_sp;

    auto ksk12_sp = ek.get<KeySwitchingKey<lvl12param>>();
    auto &ksk12 = *ksk12_sp;

    auto bkfft02_sp = ek.get<BootstrappingKeyFFT<lvl02param>>();
    auto &bkfft02 = *bkfft02_sp;

    auto bkfft12_sp = ek.get<BootstrappingKeyFFT<lvl12param>>();
    auto &bkfft12 = *bkfft12_sp;

    
    AnnihilateKey<Pout> ahk2;
    annihilatekeygen<Pout>(ahk2, sk);

    std::printf("[INFO] SecretKey & EvalKey generated.\n");

    // =====================
    
    // =====================
    autohog::SecureRng engine;
    std::uniform_int_distribution<uint32_t> binary(0, 1);

    PKE::SeededRlwePkPair<Pin>  spk_lvl1;
    PKE::SeededRlwePkPair<Pout> spk_lvl2;

    spk_lvl1 = PKE::generate_rlwe_pk_pair_seeded<Pin>(sk, engine);
    spk_lvl2 = PKE::generate_rlwe_pk_pair_seeded<Pout>(sk, engine);

    std::printf("[INFO] Seeded RLWE public key pair (lvl1 & lvl2) generated.\n");

    
    std::array<bool, Pin::n>  bits_array_lvl1{};
    std::array<bool, Pout::n> bits_array_lvl2{};

    
    TRLWE<Pin>  packed_lvl1;
    TRLWE<Pout> packed_lvl2;
    TRLWE<Pout> packed_pbs_iks;
    TRLWE<Pout> packed_pbs_only;
    TRLWE<Pout> packed_direct_iks;

    // =====================
    
    // =====================
    const int omp_threads = omp_get_max_threads();
    std::printf("\n=== Batch test NUM_SLOTS: {1,32,64,128,256,512} ===\n");
    std::printf("OpenMP threads = %d, NUM_TESTS = %d\n", omp_threads, NUM_TESTS);
    std::printf("Pin::n = %d, Pout::n = %d\n", (int)Pin::n, (int)Pout::n);

    
    std::printf("\n%-6s | %-14s | %-16s | %-14s | %-16s | %-14s | %-16s\n",
                "SLOTS",
                "avg_ms(PBS+IKS)", "mismatch(PBS+IKS)",
                "avg_ms(PBS-only)", "mismatch(PBS-only)",
                "avg_ms(DirIKS)",  "mismatch(DirIKS)");
    std::printf("---------------------------------------------------------------------------------------------------------------\n");

    for (int NUM_SLOTS : slot_list)
    {
        if (NUM_SLOTS <= 0) continue;

        
        if (NUM_SLOTS > (int)Pin::n)
        {
            std::printf("%-6d | %-14s | %-16s | %-14s | %-16s | %-14s | %-16s\n",
                        NUM_SLOTS,
                        "SKIP", "NUM_SLOTS>Pin::n",
                        "SKIP", "NUM_SLOTS>Pin::n",
                        "SKIP", "NUM_SLOTS>Pin::n");
            continue;
        }

        
        std::vector<uint8_t> m_bits(NUM_SLOTS);

        
        double total_t_pbs_iks_ms    = 0.0;
        double total_t_pbs_only_ms   = 0.0;
        double total_t_direct_iks_ms = 0.0;

        long long total_mismatch_pbs_iks  = 0;
        long long total_mismatch_pbs_only = 0;
        long long total_mismatch_direct   = 0;

        
        for (int test_id = 0; test_id < NUM_TESTS; ++test_id)
        {
            
            bits_array_lvl1.fill(false);
            bits_array_lvl2.fill(false);
            for (int i = 0; i < NUM_SLOTS; ++i)
            {
                uint8_t m = static_cast<uint8_t>(binary(engine) & 1u);
                m_bits[i] = m & 1u;
                bits_array_lvl1[i] = (m_bits[i] != 0);
                bits_array_lvl2[i] = (m_bits[i] != 0);
            }

            
            packed_lvl1 = PKE::encrypt_bits_trlwe_with_seeded_rlwe_pair<Pin>(
                spk_lvl1, bits_array_lvl1, engine);

            packed_lvl2 = PKE::encrypt_bits_trlwe_with_seeded_rlwe_pair<Pout>(
                spk_lvl2, bits_array_lvl2, engine);

            
            double t_pbs_iks_ms = HomDecompress::Decompress_PBSiks_lvl1_to_lvl2(
                packed_pbs_iks,
                packed_lvl1,
                ksk10,
                bkfft02,
                ahk2,
                NUM_SLOTS);
            total_t_pbs_iks_ms += t_pbs_iks_ms;

            {
                auto slots_pbs_iks = TRLWE_Decompress<Pout>(packed_pbs_iks, NUM_SLOTS);
                int mismatch = 0;
                for (int i = 0; i < NUM_SLOTS; ++i)
                {
                    uint8_t d = tlweSymDecrypt<Pout>(slots_pbs_iks[i], sk) & 1u;
                    if (d != m_bits[i]) mismatch++;
                }
                total_mismatch_pbs_iks += mismatch;
            }

            
            double t_pbs_only_ms = HomDecompress::Decompress_PBS_lvl1_to_lvl2(
                packed_pbs_only,
                packed_lvl1,
                bkfft12,
                ahk2,
                NUM_SLOTS);
            total_t_pbs_only_ms += t_pbs_only_ms;

            {
                auto slots_pbs_only = TRLWE_Decompress<Pout>(packed_pbs_only, NUM_SLOTS);
                int mismatch = 0;
                for (int i = 0; i < NUM_SLOTS; ++i)
                {
                    uint8_t d = tlweSymDecrypt<Pout>(slots_pbs_only[i], sk) & 1u;
                    if (d != m_bits[i]) mismatch++;
                }
                total_mismatch_pbs_only += mismatch;
            }

            
            double t_direct_iks_ms = HomDecompress::Decompress_DirectIKS_lvl1_to_lvl2(
                packed_direct_iks,
                packed_lvl1,
                ksk12,
                ahk2,
                NUM_SLOTS);
            total_t_direct_iks_ms += t_direct_iks_ms;

            {
                auto slots_direct = TRLWE_Decompress<Pout>(packed_direct_iks, NUM_SLOTS);
                int mismatch = 0;
                for (int i = 0; i < NUM_SLOTS; ++i)
                {
                    uint8_t d = tlweSymDecrypt<Pout>(slots_direct[i], sk) & 1u;
                    if (d != m_bits[i]) mismatch++;
                }
                total_mismatch_direct += mismatch;
            }
        }

        
        const double avg_t_pbs_iks_ms    = total_t_pbs_iks_ms / NUM_TESTS;
        const double avg_t_pbs_only_ms   = total_t_pbs_only_ms / NUM_TESTS;
        const double avg_t_direct_iks_ms = total_t_direct_iks_ms / NUM_TESTS;

        const long long total_bits = 1LL * NUM_SLOTS * NUM_TESTS;

        std::printf("%-6d | %-14.6f | %7lld/%-8lld | %-14.6f | %7lld/%-8lld | %-14.6f | %7lld/%-8lld\n",
                    NUM_SLOTS,
                    avg_t_pbs_iks_ms,    total_mismatch_pbs_iks,  total_bits,
                    avg_t_pbs_only_ms,   total_mismatch_pbs_only, total_bits,
                    avg_t_direct_iks_ms, total_mismatch_direct,   total_bits);
    }

    
    std::printf("\n[INFO] Done. (All times are per-call average over NUM_TESTS)\n");
    return 0;
}
