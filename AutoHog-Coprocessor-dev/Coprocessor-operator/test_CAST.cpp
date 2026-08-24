#include "HomOP.h"

#include <boost/multiprecision/cpp_int.hpp>

#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace TFHEpp;
using boost::multiprecision::cpp_int;

namespace {

static cpp_int mask_bits(int width)
{
    return (cpp_int(1) << width) - 1;
}

static std::vector<uint8_t> uint_to_bits_msb(cpp_int value, int width)
{
    value &= mask_bits(width);
    std::vector<uint8_t> bits(width);
    for (int i = 0; i < width; ++i) {
        bits[i] = static_cast<uint8_t>(((value >> (width - 1 - i)) & 1) != 0);
    }
    return bits;
}

static cpp_int bits_to_uint_msb(const std::vector<uint8_t> &bits)
{
    cpp_int value = 0;
    for (uint8_t bit : bits) {
        value <<= 1;
        value |= bit & 1u;
    }
    return value;
}

static std::vector<TLWE<lvl2param>> encrypt_uint_msb(
    cpp_int value,
    int width,
    const SecretKey &sk)
{
    const auto bits = uint_to_bits_msb(value, width);
    std::vector<TLWE<lvl2param>> out(width);
    for (int i = 0; i < width; ++i) {
        out[i] = bootsSymEncrypt<lvl2param>({bits[i]}, sk)[0];
    }
    return out;
}

static cpp_int decrypt_uint_msb(
    const std::vector<TLWE<lvl2param>> &bits,
    const SecretKey &sk)
{
    std::vector<uint8_t> plain(bits.size());
    for (size_t i = 0; i < bits.size(); ++i) {
        plain[i] = bootsSymDecrypt<lvl2param>({bits[i]}, sk)[0] & 1u;
    }
    return bits_to_uint_msb(plain);
}

static void require_eq(const char *label, const cpp_int &got, const cpp_int &expected)
{
    if (got != expected) {
        std::cout << "FAIL op=cast case=" << label
                  << " got=" << got
                  << " expected=" << expected
                  << "\n";
        throw std::runtime_error("cast test failed");
    }
}

static void require_size(const char *label, size_t got, size_t expected)
{
    if (got != expected) {
        std::cout << "FAIL op=cast case=" << label
                  << " got_size=" << got
                  << " expected_size=" << expected
                  << "\n";
        throw std::runtime_error("cast test failed");
    }
}

} // namespace

int main()
{
    try {
        SecretKey sk;
        EvalKey ek;
        ek.emplacebkfft<lvl02param>(sk);
        ek.emplaceiksk<lvl20param>(sk);
        init_constants(sk);

        const cpp_int value = cpp_int("0x8000000000000001");
        const auto v64 = encrypt_uint_msb(value, 64, sk);

        const auto zero128 = HomCast(v64, 64, 128, false);
        require_eq("zero_extend_64_to_128", decrypt_uint_msb(zero128, sk), value);

        const auto sign128 = HomCast(v64, 64, 128, true);
        const cpp_int signed_expected = ((mask_bits(64) << 64) | value) & mask_bits(128);
        require_eq("sign_extend_64_to_128", decrypt_uint_msb(sign128, sk), signed_expected);

        const auto trunc32 = HomCast(v64, 64, 32, false);
        require_eq("truncate_64_to_32", decrypt_uint_msb(trunc32, sk), value & mask_bits(32));

        const cpp_int a = cpp_int(1) << 3;
        const cpp_int b = 3;
        const auto product8 = HomMul(
            encrypt_uint_msb(a, 4, sk),
            encrypt_uint_msb(b, 4, sk),
            4,
            ek);
        require_size("mul_uses_cast_4_to_8_size", product8.size(), 8);
        require_eq("mul_uses_cast_4_to_8_value", decrypt_uint_msb(product8, sk), a * b);

        std::cout << "SUMMARY op=cast correct=5/5 semantics=zero_extend_sign_extend_truncate_mul_uses_cast\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
}
