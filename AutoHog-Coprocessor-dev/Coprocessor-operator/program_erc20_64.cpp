// Confidential 64-bit ERC-20 transfer workload.
// Circuit: LT + SUB + 2*SELECT + ADD.
// Usage: program_erc20_64 [epochs]
#include "HomOP.h"
#include "workload_common.hpp"

#include <chrono>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace TFHEpp;
using namespace autohog_workload;

int main(int argc, char **argv)
{
    try {
        constexpr int width = 64;
        const int epochs = argc > 1 ? std::stoi(argv[1]) : 1;
        if (epochs <= 0) throw std::runtime_error("epochs must be positive");
        const int threads = thread_count();

        SecretKey sk;
        EvalKey ek;
        ek.emplacebkfft<lvl02param>(sk);
        ek.emplaceiksk<lvl20param>(sk);
        init_constants(sk);

        struct Case { uint64_t sender, receiver, amount; };
        const Case cases[] = {
            {1000, 250, 300}, // accepted
            {100, 250, 300}   // rejected without revealing the branch
        };

        double total_ms = 0;
        int correct = 0;
        for (int epoch = 0; epoch < epochs; ++epoch) {
            const Case input = cases[epoch % 2];
            const auto sender = encrypt_uint(input.sender, width, sk);
            const auto receiver = encrypt_uint(input.receiver, width, sk);
            const auto amount = encrypt_uint(input.amount, width, sk);

            const auto begin = std::chrono::steady_clock::now();
            // sign=1 means insufficient funds. HomSelect(A,B,s) is s ? B : A.
            const auto sign = HomLT_unsigned_pal(sender, amount, width, ek, threads);
            const auto sender_candidate = HomSub_unsigned_pal(
                sender, amount, width, ek, false, threads).SUM;
            const auto sender_out = HomSelect(
                sender_candidate, sender, sign, width, ek, threads);
            const auto credited_amount = HomSelect(
                amount, encrypted_zero(width), sign, width, ek, threads);
            const auto receiver_out = HomAdd_unsigned_pal(
                receiver, credited_amount, width, ek, false, threads).SUM;
            const auto end = std::chrono::steady_clock::now();
            total_ms += std::chrono::duration<double, std::milli>(end - begin).count();

            const bool rejected = input.sender < input.amount;
            const uint64_t expected_sender = rejected
                ? input.sender : input.sender - input.amount;
            const uint64_t expected_receiver = rejected
                ? input.receiver : input.receiver + input.amount;
            const bool pass = decrypt_bit(sign, sk) == rejected &&
                decrypt_uint(sender_out, sk) == expected_sender &&
                decrypt_uint(receiver_out, sk) == expected_receiver;
            correct += pass;
            std::cout << "epoch=" << epoch + 1
                      << " status=" << (rejected ? "rejected" : "accepted")
                      << " sender_out=" << decrypt_uint(sender_out, sk)
                      << " receiver_out=" << decrypt_uint(receiver_out, sk)
                      << " correct=" << (pass ? "yes" : "no") << '\n';
        }

        // Each SELECT realizes sign*a + (1-sign)*b without homomorphic
        // multiplication. This is the hybrid branchless transaction circuit.
        std::cout << std::fixed << std::setprecision(3)
                  << "ERC20_64 epochs=" << epochs
                  << " avg_ms=" << total_ms / epochs
                  << " circuit=1xLT+1xSUB+2xSELECT+1xADD"
                  << " correct=" << correct << '/' << epochs << '\n';
        return correct == epochs ? 0 : 2;
    } catch (const std::exception &error) {
        std::cerr << "program_erc20_64: " << error.what() << '\n';
        return 1;
    }
}
