#include "fhe_key.h"
#include <vector>
#include <string>
#include <cstring>
#include <iostream>

using namespace std;

// --- Helper Functions ---
bool compare(const vector<uint8_t> &lhs, const vector<uint8_t>& rhs) {
    return lhs.size() == rhs.size() && memcmp(lhs.data(), rhs.data(), lhs.size()) == 0;
}

string hex_encode(const vector<uint8_t>& data) {
    string s;
    char buf[4];
    for (auto d: data) {
        sprintf(buf, "%02x", d);
        s += buf;
    }
    return s;
}

int main() {
    // 1. Initialize and generate key files (Setup Phase)
    printf("[Setup] Generating and saving keys (including AHK)...\n");
    {
        FheKey setup_key;
        setup_key.keygen(0);
        setup_key.save("client_secret.bin", "server_eval.bin", "server_ahk.bin", "server_pk.bin", "server_pk_v2.bin");
    }

    uint32_t num1 = 10;
    uint32_t num2 = 20;

    // Call to_vec from fhe_key.h directly
    vector<uint8_t> data1 = to_vec(num1);
    vector<uint8_t> data2 = to_vec(num2);

    vector<uint8_t> lhs_1, rhs_1, res_1;
    vector<uint8_t> lhs_2, rhs_2, res_2;

    // Define two independent Key objects for subsequent comparison
    FheKey key1;
    FheKey key2;

    // 3. Lambda modification: accepts an FheKey reference as a parameter
    auto compute = [&](FheKey& current_key, 
                       vector<uint8_t>& lhs_bytes, 
                       vector<uint8_t>& rhs_bytes, 
                       vector<uint8_t>& res_bytes) {
        
        // Perform load inside the Lambda to verify loading consistency
        printf("  -> Loading keys from disk (including AHK)...\n");
        current_key.load("server_eval.bin", "server_ahk.bin", "server_pk.bin", "server_pk_v2.bin");
        current_key.load_sk("client_secret.bin");

        auto lhs = current_key.trivialEncryptBatch(data1);
        auto rhs = current_key.trivialEncryptBatch(data2);

        auto lhs_ct = decompress_ciphertexts(lhs);
        auto rhs_ct = decompress_ciphertexts(rhs);
        
        // Perform homomorphic addition (Consumes RNG state)
        auto res_ct = current_key.add(lhs_ct, rhs_ct);

        lhs_bytes = serialize_compressed_ciphertexts(lhs);
        rhs_bytes = serialize_compressed_ciphertexts(rhs);
        res_bytes = serialize_ciphertexts(res_ct);

        auto pt = current_key.decryptBatch(res_ct);
        printf("  -> Decrypted Result: %s\n", hex_encode(pt).c_str());
    };

    // --- Round 1 (Using Key Instance 1) ---
    printf("\n--- Round 1 (Using Key Instance 1) ---\n");
    compute(key1, lhs_1, rhs_1, res_1);

    // --- Round 2 (Using Key Instance 2) ---
    printf("\n--- Round 2 (Using Key Instance 2) ---\n");
    compute(key2, lhs_2, rhs_2, res_2);

    // ==========================================
    // 4. Core Comparison Logic
    // ==========================================
    printf("\n=== Consistency Report ===\n");

    // 4.1 Compare if two Key objects are equal (calls operator==)
    // This triggers the serialized comparison defined in fhe_key.h
    printf("Checking Key Equality (Key1 vs Key2)...\n");
    if (key1 == key2) {
        printf("SUCCESS Key Objects are IDENTICAL. (Load is deterministic)\n");
    } else {
        printf("FAILED Key Objects DIFFER. (Load introduced randomness or metadata change)\n");
    }

    // 4.2 Compare calculation results
    bool lhs_cmp = compare(lhs_1, lhs_2);
    bool rhs_cmp = compare(rhs_1, rhs_2);
    bool res_cmp = compare(res_1, res_2);
    
    printf("\nChecking Ciphertext Output Consistency:\n");
    printf("LHS Identical? : %d\n", lhs_cmp);
    printf("RHS Identical? : %d\n", rhs_cmp);
    printf("RES Identical? : %d\n", res_cmp);
    
    if (res_cmp && (key1 == key2)) {
        printf("\n SUCCESS Conclusion: System has perfect determinism (Keys match AND Results match).\n");
    } else {
        printf("\n FAILED Conclusion: Inconsistency detected.\n");
    }

    return 0;
}
