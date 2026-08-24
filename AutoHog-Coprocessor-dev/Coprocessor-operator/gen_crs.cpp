#include <array>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include <openssl/rand.h>

#include "ZKProof.h"

namespace {

void secure_zero(void* ptr, size_t len)
{
    volatile unsigned char* p = static_cast<volatile unsigned char*>(ptr);
    while (len-- > 0) {
        *p++ = 0;
    }
}

class SeedBuffer {
public:
    std::array<uint8_t, ZKProof::SEED_LEN>& get() { return bytes_; }
    const std::array<uint8_t, ZKProof::SEED_LEN>& get() const { return bytes_; }

    ~SeedBuffer() { secure_zero(bytes_.data(), bytes_.size()); }

private:
    std::array<uint8_t, ZKProof::SEED_LEN> bytes_{};
};

void fill_seed(std::array<uint8_t, ZKProof::SEED_LEN>& seed)
{
    if (RAND_bytes(seed.data(), static_cast<int>(seed.size())) != 1) {
        throw std::runtime_error("RAND_bytes failed while generating CRS seed");
    }
}

void print_usage(const char* argv0)
{
    std::cout
        << "Usage: " << argv0
        << " [output_path] [--force] --allow-insecure-single-party\n"
        << "  output_path: CRS file path (default: crs_lvl2.bin)\n"
        << "  --force:     overwrite an existing CRS file\n"
        << "  --allow-insecure-single-party:\n"
        << "               acknowledge that this test/development generator creates\n"
        << "               toxic waste in one process and is not a production ceremony\n";
}

} // namespace

int main(int argc, char** argv) try
{
    using P = TFHEpp::lvl2param;

    std::string output_path = "crs_lvl2.bin";
    std::string hash_output_path;
    bool force_overwrite = false;
    bool allow_insecure_single_party = false;
    bool saw_output_path = false;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            print_usage(argv[0]);
            return 0;
        }
        if (arg == "--force") {
            force_overwrite = true;
            continue;
        }
        if (arg == "--allow-insecure-single-party") {
            allow_insecure_single_party = true;
            continue;
        }
        if (saw_output_path) {
            std::cerr << "Unexpected argument: " << arg << "\n";
            print_usage(argv[0]);
            return 1;
        }
        output_path = arg;
        saw_output_path = true;
    }

    hash_output_path = output_path + ".sha256";

    if (!allow_insecure_single_party) {
        std::cerr
            << "Refusing single-party CRS generation by default.\n"
            << "Production must use an audited multi-party ceremony and a pinned artifact digest.\n"
            << "For tests only, re-run with --allow-insecure-single-party.\n";
        return 2;
    }

    if (std::filesystem::exists(output_path) && !force_overwrite) {
        std::cerr
            << "CRS file already exists: " << output_path << "\n"
            << "Re-run with --force to overwrite it.\n";
        return 1;
    }

    constexpr uint64_t ZK_B = 1ULL << 16;
    constexpr uint64_t ZK_BR = static_cast<uint64_t>(P::n / 2 + 1);
    constexpr uint64_t ZK_T = 32;
    constexpr uint64_t ZK_MSBS = 0;
    constexpr size_t ZK_EFFECTIVE_T = ZK_T >> ZK_MSBS;
    constexpr size_t ZK_STEP_T = 5;
    constexpr size_t ZK_STEP_B = 17;
    constexpr size_t ZK_STEP_BR = 11;
    static_assert(ZK_EFFECTIVE_T == 32);

    std::cout << "=== INSECURE SINGLE-PARTY CRS GENERATION (TEST ONLY) ===\n";
    std::cout << "WARNING: do not deploy this CRS in production.\n";
    std::cout << "Output:              " << output_path << "\n";
    std::cout << "Parameter set:       lvl2param\n";
    std::cout << "RLWE dimension d:    " << P::n << "\n";
    std::cout << "Statement dimension: " << P::n << "\n";
    std::cout << "Bounds:              B=" << ZK_B << ", Br=" << ZK_BR << "\n";
    std::cout << "Encoding:            q=2^64, t=" << ZK_T << ", msbs=" << ZK_MSBS << "\n\n";

    SeedBuffer alpha_seed;
    fill_seed(alpha_seed.get());

    const auto start = std::chrono::high_resolution_clock::now();
    auto params = ZKProof::generate_public_params<P>(
        0,
        ZK_B,
        ZK_BR,
        0,
        ZK_T,
        ZK_MSBS,
        alpha_seed.get()
    );
    const size_t ZK_STATEMENT_K = params.k;
    const size_t ZK_EXPECTED_N =
        P::n + ZK_STATEMENT_K * ZK_STEP_T + (P::n + ZK_STATEMENT_K) * (ZK_STEP_B + ZK_STEP_BR) + 1;
    ZKProof::save_public_params<P>(params, output_path);
    const std::string digest = ZKProof::compute_file_sha256_hex(output_path);
    {
        std::ofstream ofs(hash_output_path, std::ios::binary);
        if (!ofs) {
            throw std::runtime_error("Cannot open hash file for writing: " + hash_output_path);
        }
        ofs << digest << '\n';
        if (!ofs) {
            throw std::runtime_error("Failed writing hash file: " + hash_output_path);
        }
    }
    auto loaded_params = ZKProof::load_public_params<P>(output_path);
    const auto end = std::chrono::high_resolution_clock::now();

    if (loaded_params.d != P::n ||
        loaded_params.k != ZK_STATEMENT_K ||
        loaded_params.n != ZK_EXPECTED_N ||
        loaded_params.b != ZK_B ||
        loaded_params.b_r != ZK_BR ||
        loaded_params.q != 0 ||
        loaded_params.t != ZK_T ||
        loaded_params.msbs_zero_padding_bit_count != ZK_MSBS) {
        throw std::runtime_error("Saved CRS failed round-trip validation");
    }

    const auto elapsed_ms =
        std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(end - start);
    const auto file_size = std::filesystem::file_size(output_path);

    std::cout << "CRS generated and reloaded successfully.\n";
    std::cout << "Elapsed: " << elapsed_ms.count() << " ms\n";
    std::cout << "File size: " << file_size << " bytes\n";
    std::cout << "SHA-256: " << digest << "\n";
    std::cout << "SHA-256 file: " << hash_output_path << "\n";
    std::cout << "Use the pinned digest via --crs-sha256 or AUTOHOG_CRS_SHA256.\n";
    std::cout << "Seed and derived scalar buffers are wiped after use.\n";
    std::cout << "This memory hygiene does not replace a multi-party ceremony.\n";

    return 0;
}
catch (const std::exception& e)
{
    std::cerr << "gen_crs failed: " << e.what() << "\n";
    return 1;
}
