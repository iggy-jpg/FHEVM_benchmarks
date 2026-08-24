#pragma once

#include "Autohog_Cop.hpp"

#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cstdint>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace autohog_workload {

using TFHEpp::SecretKey;
using TFHEpp::TLWE;
using TFHEpp::lvl2param;
using boost::multiprecision::cpp_int;

inline int thread_count()
{
#ifdef _OPENMP
    return std::max(1, omp_get_max_threads());
#else
    return 1;
#endif
}

inline cpp_int bit_mask(int width) { return (cpp_int(1) << width) - 1; }

inline std::vector<uint8_t> to_bits(cpp_int value, int width)
{
    value &= bit_mask(width);
    std::vector<uint8_t> bits(width);
    for (int i = 0; i < width; ++i)
        bits[i] = static_cast<uint8_t>(((value >> (width - 1 - i)) & 1) != 0);
    return bits;
}

inline std::vector<TLWE<lvl2param>> encrypt_uint(
    cpp_int value, int width, const SecretKey &sk)
{
    const auto bits = to_bits(value, width);
    std::vector<TLWE<lvl2param>> out(width);
    for (int i = 0; i < width; ++i)
        out[i] = TFHEpp::bootsSymEncrypt<lvl2param>({bits[i]}, sk)[0];
    return out;
}

inline cpp_int decrypt_uint(
    const std::vector<TLWE<lvl2param>> &value, const SecretKey &sk)
{
    cpp_int out = 0;
    for (const auto &bit : value)
        out = (out << 1) |
              (TFHEpp::bootsSymDecrypt<lvl2param>({bit}, sk)[0] & 1);
    return out;
}

inline int decrypt_bit(const TLWE<lvl2param> &value, const SecretKey &sk)
{
    return TFHEpp::bootsSymDecrypt<lvl2param>({value}, sk)[0] & 1;
}

inline std::vector<TLWE<lvl2param>> encrypted_zero(int width)
{
    return std::vector<TLWE<lvl2param>>(width, TFHEpp::CT_ZERO);
}

} // namespace autohog_workload
