#include <cassert>
#include <chrono> 
#include <iostream>
#include <random>
#include <vector>            
#include <tfhe++.hpp>

int main()
{
    constexpr uint32_t num_test = 1000;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> binary(0, 1);

    TFHEpp::SecretKey sk;
    TFHEpp::EvalKey ek;
    ek.emplacebkfft<TFHEpp::lvl02param>(sk);
    ek.emplaceiksk<TFHEpp::lvl20param>(sk);

    std::chrono::duration<double, std::milli> total_time(0);

    for (uint32_t test = 0; test < num_test; test++) {
        uint8_t p = static_cast<uint8_t>(binary(engine) > 0);

        std::vector<uint8_t> plain_vec{p};
        auto enc_vec = TFHEpp::bootsSymEncrypt<TFHEpp::lvl2param>(plain_vec, sk);
        TFHEpp::TLWE<TFHEpp::lvl2param> tlwe = enc_vec.at(0);

        TFHEpp::TLWE<TFHEpp::lvl2param> bootedtlwe;

        auto start = std::chrono::high_resolution_clock::now();
        TFHEpp::GateBootstrapping<
            TFHEpp::lvl20param,
            TFHEpp::lvl02param,
            TFHEpp::lvl2param::μ
        >(bootedtlwe, tlwe, ek);
        auto end = std::chrono::high_resolution_clock::now();

        auto diff = std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(end - start);
        total_time += diff;

        auto dec_before_vec = TFHEpp::bootsSymDecrypt<TFHEpp::lvl2param>(
            std::vector<TFHEpp::TLWE<TFHEpp::lvl2param>>{tlwe}, sk);
        auto dec_after_vec = TFHEpp::bootsSymDecrypt<TFHEpp::lvl2param>(
            std::vector<TFHEpp::TLWE<TFHEpp::lvl2param>>{bootedtlwe}, sk);

        uint8_t p_before = dec_before_vec.at(0);
        uint8_t p_after  = dec_after_vec.at(0);

        assert(p == p_before);
        assert(p == p_after);
    }

    std::cout << "Total bootstrapping time: " << total_time.count() << " ms" << std::endl;
    std::cout << "Average per operation: " << total_time.count() / num_test << " ms" << std::endl;

    std::cout << "Passed" << std::endl;
    return 0;
}
