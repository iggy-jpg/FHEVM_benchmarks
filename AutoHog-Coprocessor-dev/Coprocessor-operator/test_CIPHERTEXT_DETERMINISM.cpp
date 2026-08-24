#include "HomOP.h"

#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
#include <openssl/evp.h>

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace TFHEpp;

namespace {

constexpr std::uint32_t kFormatVersion = 1;

struct InputArtifact {
    std::uint32_t version = kFormatVersion;
    std::uint32_t width = 8;
    std::uint64_t clear_a = 0x5a;
    std::uint64_t clear_b = 0x13;
    SecretKey secret_key;
    EvalKey coefficient_eval_key;
    TLWE<lvl2param> ct_zero{};
    TLWE<lvl2param> ct_one{};
    std::vector<TLWE<lvl2param>> encrypted_a;
    std::vector<TLWE<lvl2param>> encrypted_b;

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(version, width, clear_a, clear_b, secret_key, coefficient_eval_key,
                ct_zero, ct_one, encrypted_a, encrypted_b);
    }
};

struct OutputArtifact {
    std::uint32_t version = kFormatVersion;
    std::uint32_t width = 0;
    std::string backend;
    int threads = 1;
    bool semantic_ok = false;
    std::vector<TLWE<lvl2param>> add;
    std::vector<TLWE<lvl2param>> gt;
    std::vector<TLWE<lvl2param>> mul;
    std::vector<TLWE<lvl2param>> quotient;
    std::vector<TLWE<lvl2param>> remainder;

    template <class Archive>
    void serialize(Archive &archive)
    {
        archive(version, width, backend, threads, semantic_ok, add, gt, mul,
                quotient, remainder);
    }
};

template <class T>
void save_portable(const T &value, const std::string &path)
{
    std::ofstream stream(path, std::ios::binary);
    if (!stream) throw std::runtime_error("cannot write " + path);
    cereal::PortableBinaryOutputArchive archive(stream);
    archive(value);
}

template <class T>
T load_portable(const std::string &path)
{
    std::ifstream stream(path, std::ios::binary);
    if (!stream) throw std::runtime_error("cannot read " + path);
    cereal::PortableBinaryInputArchive archive(stream);
    T value;
    archive(value);
    return value;
}

std::string sha256_file(const std::string &path)
{
    std::ifstream stream(path, std::ios::binary);
    if (!stream) throw std::runtime_error("cannot hash " + path);
    EVP_MD_CTX *ctx = EVP_MD_CTX_new();
    if (ctx == nullptr) throw std::runtime_error("EVP_MD_CTX_new failed");
    unsigned char digest[EVP_MAX_MD_SIZE];
    unsigned int digest_size = 0;
    std::vector<char> buffer(1 << 20);
    bool ok = EVP_DigestInit_ex(ctx, EVP_sha256(), nullptr) == 1;
    while (ok && stream) {
        stream.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
        const auto count = stream.gcount();
        if (count > 0)
            ok = EVP_DigestUpdate(ctx, buffer.data(), static_cast<std::size_t>(count)) == 1;
    }
    ok = ok && EVP_DigestFinal_ex(ctx, digest, &digest_size) == 1;
    EVP_MD_CTX_free(ctx);
    if (!ok) throw std::runtime_error("SHA-256 failed");
    std::ostringstream out;
    for (unsigned int i = 0; i < digest_size; ++i)
        out << std::hex << std::setw(2) << std::setfill('0') << unsigned(digest[i]);
    return out.str();
}

std::vector<TLWE<lvl2param>> encrypt_value(std::uint64_t value, int width,
                                            const SecretKey &sk)
{
    std::vector<TLWE<lvl2param>> result(width);
    for (int i = 0; i < width; ++i) {
        const std::uint8_t bit = (value >> (width - 1 - i)) & 1U;
        result[i] = bootsSymEncrypt<lvl2param>({bit}, sk)[0];
    }
    return result;
}

std::uint64_t decrypt_value(const std::vector<TLWE<lvl2param>> &value,
                            const SecretKey &sk)
{
    std::uint64_t result = 0;
    for (const auto &bit : value)
        result = (result << 1) | (bootsSymDecrypt<lvl2param>({bit}, sk)[0] & 1U);
    return result;
}

std::string backend_name()
{
#if defined(USE_FFTW3)
    return "fftw3";
#elif defined(USE_MKL)
    return "mkl";
#elif defined(USE_SPQLIOX_AARCH64)
    return "spqliox_aarch64";
#elif defined(USE_CONCRETE_FFT)
    return "concrete_fft";
#elif defined(USE_SPQLIOS_ARITHMETIC)
    return "spqlios_arithmetic";
#else
    return "spqlios";
#endif
}

struct Difference {
    std::uint64_t coefficients = 0;
    std::uint64_t sanitized4 = 0;
    std::size_t first_ciphertext = 0;
    std::size_t first_coefficient = 0;
    int lowest_differing_bit = 64;
    int highest_differing_bit = -1;
};

int lowest_bit(std::uint64_t value)
{
    int bit = 0;
    while ((value & 1U) == 0U) { value >>= 1; ++bit; }
    return bit;
}

int highest_bit(std::uint64_t value)
{
    int bit = 0;
    while (value >>= 1) ++bit;
    return bit;
}

Difference compare_vectors(const std::vector<TLWE<lvl2param>> &lhs,
                           const std::vector<TLWE<lvl2param>> &rhs)
{
    if (lhs.size() != rhs.size()) throw std::runtime_error("ciphertext vector size mismatch");
    Difference difference;
    bool first = true;
    for (std::size_t c = 0; c < lhs.size(); ++c) {
        for (std::size_t i = 0; i < lhs[c].size(); ++i) {
            const std::uint64_t xor_value = lhs[c][i] ^ rhs[c][i];
            if (xor_value == 0) continue;
            if (first) {
                difference.first_ciphertext = c;
                difference.first_coefficient = i;
                first = false;
            }
            ++difference.coefficients;
            if ((lhs[c][i] & ~std::uint64_t{0xf}) !=
                (rhs[c][i] & ~std::uint64_t{0xf}))
                ++difference.sanitized4;
            difference.lowest_differing_bit =
                std::min(difference.lowest_differing_bit, lowest_bit(xor_value));
            difference.highest_differing_bit =
                std::max(difference.highest_differing_bit, highest_bit(xor_value));
        }
    }
    return difference;
}

void print_difference(const std::string &name, const Difference &difference)
{
    std::cout << "[COMPARE] op=" << name
              << " coefficient_mismatches=" << difference.coefficients
              << " sanitized4_mismatches=" << difference.sanitized4;
    if (difference.coefficients != 0)
        std::cout << " first_ciphertext=" << difference.first_ciphertext
                  << " first_coefficient=" << difference.first_coefficient
                  << " lowest_differing_bit=" << difference.lowest_differing_bit
                  << " highest_differing_bit=" << difference.highest_differing_bit;
    std::cout << '\n';
}

} // namespace

int main(int argc, char **argv)
{
    try {
        if (argc < 3) {
            std::cerr << "usage: test_CIPHERTEXT_DETERMINISM "
                         "generate <artifact> [width] | evaluate <artifact> <output> [threads] | "
                         "compare <output-a> <output-b>\n";
            return 2;
        }
        const std::string mode = argv[1];
        if (mode == "generate") {
            const std::string path = argv[2];
            const int width = argc > 3 ? std::stoi(argv[3]) : 8;
            if (width < 4 || width > 32) throw std::runtime_error("width must be in [4, 32]");
            InputArtifact artifact;
            artifact.width = width;
            artifact.clear_a &= (std::uint64_t{1} << width) - 1;
            artifact.clear_b &= (std::uint64_t{1} << width) - 1;
            artifact.coefficient_eval_key.emplacebk<lvl02param>(artifact.secret_key);
            artifact.coefficient_eval_key.emplaceiksk<lvl20param>(artifact.secret_key);
            init_constants(artifact.secret_key);
            artifact.ct_zero = CT_ZERO;
            artifact.ct_one = CT_ONE;
            artifact.encrypted_a = encrypt_value(artifact.clear_a, width, artifact.secret_key);
            artifact.encrypted_b = encrypt_value(artifact.clear_b, width, artifact.secret_key);
            const auto duplicate = encrypt_value(artifact.clear_a, width, artifact.secret_key);
            const auto randomization = compare_vectors(artifact.encrypted_a, duplicate);
            save_portable(artifact, path);
            std::cout << "[GENERATE] artifact=" << path << " sha256=" << sha256_file(path)
                      << " same_key_reencrypt_mismatches=" << randomization.coefficients << '\n';
            return randomization.coefficients == 0 ? 1 : 0;
        }
        if (mode == "evaluate") {
            if (argc < 4) throw std::runtime_error("evaluate requires artifact and output");
            const int threads = argc > 4 ? std::stoi(argv[4]) : 1;
#ifdef _OPENMP
            omp_set_dynamic(0);
            omp_set_num_threads(threads);
#endif
            InputArtifact input = load_portable<InputArtifact>(argv[2]);
            if (input.version != kFormatVersion) throw std::runtime_error("artifact version mismatch");
            input.coefficient_eval_key.emplacebk2bkfft<lvl02param>();
            CT_ZERO = input.ct_zero;
            CT_ONE = input.ct_one;

            OutputArtifact output;
            output.width = input.width;
            output.backend = backend_name();
            output.threads = threads;
            output.add = HomAdd_pal(input.encrypted_a, input.encrypted_b, input.width,
                                    input.coefficient_eval_key, false, threads).SUM;
            output.gt = {HomGT_pal(input.encrypted_a, input.encrypted_b, input.width,
                                   input.coefficient_eval_key, threads)};
            output.mul = HomMul_pal(input.encrypted_a, input.encrypted_b, input.width,
                                    input.coefficient_eval_key, threads);
            auto division = HomDiv_pal(input.encrypted_a, input.encrypted_b, input.width,
                                       input.coefficient_eval_key, threads);
            output.quotient = std::move(division.quotient);
            output.remainder = std::move(division.remainder);
            const std::uint64_t mask = (std::uint64_t{1} << input.width) - 1;
            const auto signed_value = [width = input.width](std::uint64_t value) {
                const std::uint64_t sign = std::uint64_t{1} << (width - 1);
                return static_cast<std::int64_t>((value ^ sign) - sign);
            };
            output.semantic_ok =
                decrypt_value(output.add, input.secret_key) ==
                    ((input.clear_a + input.clear_b) & mask) &&
                decrypt_value(output.gt, input.secret_key) ==
                    (signed_value(input.clear_a) > signed_value(input.clear_b)) &&
                decrypt_value(output.mul, input.secret_key) == input.clear_a * input.clear_b &&
                decrypt_value(output.quotient, input.secret_key) == input.clear_a / input.clear_b &&
                decrypt_value(output.remainder, input.secret_key) == input.clear_a % input.clear_b;
            save_portable(output, argv[3]);
            std::cout << "[EVALUATE] backend=" << output.backend << " threads=" << threads
                      << " semantic_ok=" << output.semantic_ok << " output=" << argv[3]
                      << " sha256=" << sha256_file(argv[3]) << '\n';
            return output.semantic_ok ? 0 : 1;
        }
        if (mode == "compare") {
            if (argc < 4) throw std::runtime_error("compare requires two outputs");
            const auto lhs = load_portable<OutputArtifact>(argv[2]);
            const auto rhs = load_portable<OutputArtifact>(argv[3]);
            if (lhs.width != rhs.width) throw std::runtime_error("output width mismatch");
            std::uint64_t raw_mismatches = 0;
            for (const auto &[name, difference] : std::vector<std::pair<std::string, Difference>>{
                     {"add", compare_vectors(lhs.add, rhs.add)},
                     {"gt", compare_vectors(lhs.gt, rhs.gt)},
                     {"mul", compare_vectors(lhs.mul, rhs.mul)},
                     {"div_quotient", compare_vectors(lhs.quotient, rhs.quotient)},
                     {"div_remainder", compare_vectors(lhs.remainder, rhs.remainder)}}) {
                print_difference(name, difference);
                raw_mismatches += difference.coefficients;
            }
            std::cout << "[COMPARE_RESULT] lhs_backend=" << lhs.backend
                      << " rhs_backend=" << rhs.backend
                      << " semantic_ok=" << (lhs.semantic_ok && rhs.semantic_ok)
                      << " byte_exact=" << (raw_mismatches == 0)
                      << " coefficient_mismatches=" << raw_mismatches << '\n';
            return raw_mismatches == 0 && lhs.semantic_ok && rhs.semantic_ok ? 0 : 1;
        }
        throw std::runtime_error("unknown mode: " + mode);
    } catch (const std::exception &error) {
        std::cerr << "[FAIL] " << error.what() << '\n';
        return 1;
    }
}
