#include "fhe_key.h"
#include <vector>
#include <iostream>
#include <iomanip>
#include <cstring>

using namespace std;
using namespace TFHEpp;

// ==========================================
// 0. Helper Functions
// ==========================================

// Compare if two vectors are completely identical
bool compare(const vector<uint8_t> &lhs, const vector<uint8_t>& rhs) {
    return lhs.size() == rhs.size() && memcmp(lhs.data(), rhs.data(), lhs.size()) == 0;
}

// Convert vector to Hex string
string hex_encode(const vector<uint8_t>& data) {
    string s;
    char buf[4];
    for (auto d: data) {
        sprintf(buf, "%02x", d);
        s += buf;
    }
    return s;
}

// Integer to Vector (Big Endian) - Ensure deterministic conversion
Bytes to_vec_u32(uint32_t val) {
    Bytes result(4);
    result[0] = (val >> 24) & 0xFF;
    result[1] = (val >> 16) & 0xFF;
    result[2] = (val >> 8) & 0xFF;
    result[3] = val & 0xFF;
    return result;
}

// Check if the ciphertext contains randomness (Check if Mask is non-zero)
bool check_confidentiality(const tlwe_vec_ptr& ct) {
    // Trivial ciphertext has all-zero Masks.
    // Standard ciphertext should have random Masks.
    for (const auto& sample : ct->data) {
        for (int i = 0; i < lvl2param::n; ++i) { 
            if (sample[i] != 0) return true; // Found non-zero mask, implies noise exists
        }
    }
    return false; // All masks are zero, implies naked Trivial ciphertext
}

// ==========================================
// 1. Main Test Flow
// ==========================================

int main(int argc, const char* argv[]) {
    // 0. Initialization
    printf("[Init] Generating Keys...\n");
    {
        FheKey key;
        key.keygen(0);
        key.save("client_secret.bin", "server_eval.bin", "server_pk.bin", "server_pk_v2.bin");
    }

    // Prepare data: 30 - 10 = 20
    uint32_t num1 = 30; // 00011110
    uint32_t num2 = 10; // 00001010
    vector<uint8_t> data1 = to_vec_u32(num1);
    vector<uint8_t> data2 = to_vec_u32(num2);

    // Store results from two runs for comparison
    vector<uint8_t> bin_sub_1, bin_sub_2;      
    vector<uint8_t> bin_std_1, bin_std_2;      
    vector<uint8_t> bin_rec_1, bin_rec_2;

    auto compute = [&](int run_id, 
                       vector<uint8_t>& out_sub, 
                       vector<uint8_t>& out_std,
                       vector<uint8_t>& out_rec) {
        printf("\n--- Run ID: %d ---\n", run_id);
        
        FheKey key;
        // Load keys (Assuming annihilation key generation from SK if file not present)
        key.load("server_eval.bin", "server_ahk.bin", "server_pk.bin", "server_pk_v2.bin");
        key.load_sk("client_secret.bin");

        // Prepare Trivial Inputs
        auto lhs = key.trivialEncryptBatch(data1);
        auto rhs = key.trivialEncryptBatch(data2);
        auto lhs_ct = decompress_ciphertexts(lhs);
        auto rhs_ct = decompress_ciphertexts(rhs);

        // ------------------------------------------------
        // [Test 1] Deterministic Subtraction
        // ------------------------------------------------
        // Use member function key.trivialSub
        auto res_sub = key.trivialSub(lhs_ct, rhs_ct);
        
        auto pt_sub = key.decryptBatch(res_sub);
        printf("[Test 1] Sub (30 - 10): %s (Expect 14)\n", hex_encode(pt_sub).c_str());
        out_sub = serialize_ciphertexts(res_sub);

        // ------------------------------------------------
        // [Test 2] Trivial -> Standard (Inject Noise)
        // ------------------------------------------------
        // Use member function key.trivialToStandard
        auto res_std = key.trivialToStandard(res_sub);
        
        auto pt_std = key.decryptBatch(res_std);
        bool is_secure = check_confidentiality(res_std);
        
        printf("[Test 2] Triv->Std Decrypt: %s (Expect 14)\n", hex_encode(pt_std).c_str());
        printf("[Test 2] Confidentiality:   %s\n", is_secure ? "PASS (Has Noise)" : "FAIL (No Noise)");
        
        out_std = serialize_ciphertexts(res_std);

        // ------------------------------------------------
        // [Test 3] Standard -> Trivial (Clean Consensus)
        // ------------------------------------------------
        // Use member function key.standardToTrivial
        auto res_rec = key.standardToTrivial(res_std);

        auto pt_rec = key.decryptBatch(res_rec);
        printf("[Test 3] Std->Triv Decrypt: %s (Expect 14)\n", hex_encode(pt_rec).c_str());

        out_rec = serialize_ciphertexts(res_rec);
    };

    // Run twice
    compute(1, bin_sub_1, bin_std_1, bin_rec_1);
    compute(2, bin_sub_2, bin_std_2, bin_rec_2);

    // ==========================================
    // Final Consistency Report
    // ==========================================
    printf("\n=== Final Consistency Report ===\n");

    // 1. Trivial Subtraction (Must be Identical)
    bool sub_ok = compare(bin_sub_1, bin_sub_2);
    printf("1. Trivial Sub Determinism:      %s %s\n", 
           sub_ok ? "YES" : "NO", sub_ok ? "SUCCESS" : "FAIL");

    // 2. Standard Ciphertext (Must differ due to randomness)
    bool std_ok = !compare(bin_std_1, bin_std_2); 
    // Note: compare returns true if identical, so we negate it.
    printf("2. Triv->Std Randomness:         %s %s\n", 
           std_ok ? "YES" : "NO", std_ok ? "SUCCESS (Secure)" : "FAIL (Insecure: Deterministic)");

    // 3. Recovered Trivial (Must be Identical, Consensus Reached)
    bool rec_ok = compare(bin_rec_1, bin_rec_2);
    printf("3. Std->Triv Consensus Cleaning: %s %s\n", 
           rec_ok ? "YES" : "NO", rec_ok ? "SUCCESS" : "FAIL");

    return 0;
}
