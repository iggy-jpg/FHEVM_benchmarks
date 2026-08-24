#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <vector>                           
#include "HomPKE.h"
#include "SaveFile.h"
#include "SecureRandom.h"

int main()
{
    constexpr uint32_t num_test = 1;

    autohog::SecureRng engine;
    std::uniform_int_distribution<uint32_t> binary(0, 1);

    TFHEpp::SecretKey sk;
    TFHEpp::EvalKey ek;
    ek.emplacebkfft<TFHEpp::lvl02param>(sk);
    ek.emplaceiksk<TFHEpp::lvl20param>(sk);
    // auto pk = PKE::generate_public_key<TFHEpp::lvl2param>(sk);
    auto cpk = PKE::generate_compressed_public_key<TFHEpp::lvl2param>(sk);

    // TFHEpp::savePublicKey<TFHEpp::lvl2param>(pk, "publickey.bin");
    // // TFHEpp::printFileSize("publickey.bin");
    // TFHEpp::loadPublicKey<TFHEpp::lvl2param>(pk, "publickey.bin");

    TFHEpp::saveCompressedPublicKey<TFHEpp::lvl2param>(cpk, "publickey.bin");
    TFHEpp::printFileSize("publickey.bin");
    TFHEpp::loadCompressedPublicKey<TFHEpp::lvl2param>(cpk, "publickey.bin");

    std::chrono::duration<double, std::milli> total_boot_time(0.0);
    std::chrono::duration<double, std::milli> total_pke_time(0.0);

    for (uint32_t t = 0; t < num_test; ++t) {
        uint8_t m = static_cast<uint8_t>(binary(engine) > 0);

        auto pke_start = std::chrono::high_resolution_clock::now();
        // TFHEpp::TLWE<TFHEpp::lvl2param> ct_pke =
        //     PKE::encrypt_public_key_bit<TFHEpp::lvl2param>(pk, m, engine);、
        TFHEpp::TLWE<TFHEpp::lvl2param> ct_pke =
               PKE::encrypt_with_compressed_public_key_bit<TFHEpp::lvl2param>(cpk, m, engine);
        auto pke_end = std::chrono::high_resolution_clock::now();
        total_pke_time += std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(pke_end - pke_start);

        TFHEpp::TLWE<TFHEpp::lvl2param> booted;
        auto boot_start = std::chrono::high_resolution_clock::now();

        
        TFHEpp::GateBootstrapping<
            TFHEpp::lvl20param,        
            TFHEpp::lvl02param,        
            TFHEpp::lvl2param::μ       
        >(booted, ct_pke, ek);          // [CHANGED]

        auto boot_end = std::chrono::high_resolution_clock::now();
        total_boot_time += std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(boot_end - boot_start);

        
        uint8_t dec1 = TFHEpp::bootsSymDecrypt<TFHEpp::lvl2param>(
            std::vector<TFHEpp::TLWE<TFHEpp::lvl2param>>{ ct_pke },  
            sk
        ).at(0);                                                    

        
        uint8_t dec2 = TFHEpp::bootsSymDecrypt<TFHEpp::lvl2param>(
            std::vector<TFHEpp::TLWE<TFHEpp::lvl2param>>{ booted },
            sk
        ).at(0);

        std::cout << "m=" << int(m) << " dec(PKE)=" << int(dec1) << " dec(PBS)=" << int(dec2) << "\n";
        assert(m == dec1);
        assert(m == dec2);
    }

    // std::cout << "Total PKE time: " << total_pke_time.count() << " ms\n";
    // std::cout << "Avg PKE per op: " << (total_pke_time.count() / num_test) << " ms\n";
    // std::cout << "Total bootstrapping time: " << total_boot_time.count() << " ms\n";
    // std::cout << "Avg bootstrapping per op: " << (total_boot_time.count() / num_test) << " ms\n";
    // std::cout << "Passed\n";
    return 0;
}
