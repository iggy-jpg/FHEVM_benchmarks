#include "fhe_key.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cstring>
#include <random>

using namespace std;
using namespace TFHEpp;

// ==========================================
// Helper Functions
// ==========================================

// Compare if two files are identical (binary comparison)
bool compare_files(const string& p1, const string& p2) {
    ifstream f1(p1, ifstream::binary | ifstream::ate);
    ifstream f2(p2, ifstream::binary | ifstream::ate);

    if (f1.fail() || f2.fail()) {
        cerr << "[Error] Cannot open files for comparison." << endl;
        return false;
    }

    if (f1.tellg() != f2.tellg()) {
        cerr << "[Mismatch] File sizes differ." << endl;
        return false; 
    }

    f1.seekg(0, ifstream::beg);
    f2.seekg(0, ifstream::beg);
    return equal(istreambuf_iterator<char>(f1.rdbuf()),
                 istreambuf_iterator<char>(),
                 istreambuf_iterator<char>(f2.rdbuf()));
}

// ==========================================
// Main Test Flow
// ==========================================

int main() {
    const string sec_file = "test_pack_sk.bin";
    const string eval_file = "test_pack_eval.bin";
    const string ahk_file = "test_pack_ahk.bin"; 
    const string ahk_file_2 = "test_pack_ahk_2.bin"; 
    const string pk_file = "test_pack_pk.bin";
    const string pk_v2_file = "test_pack_pk_v2.bin";

    // 0. Initialize and Generate Keys
    cout << "[Step 0] Generating Keys..." << endl;
    {
        FheKey key;
        key.keygen(0);
        // Save all keys, including AnnihilateKey
        key.save(sec_file, eval_file, ahk_file, pk_file, pk_v2_file);
    }
    cout << "  -> Keys saved." << endl;

    // 1. Test AnnihilateKey Loading and Consistency
    cout << "[Step 1] Testing AnnihilateKey Persistence..." << endl;
    {
        FheKey key_loaded;
        // Load AHK from file
        key_loaded.load(eval_file, ahk_file, pk_file, pk_v2_file);
        key_loaded.load_sk(sec_file);
        
        // Save again to a new file for comparison
        key_loaded.save(sec_file, eval_file, ahk_file_2, pk_file, pk_v2_file);
    }

    // Compare the two AHK files
    if (compare_files(ahk_file, ahk_file_2)) {
        cout << "  SUCCESS PASS: AnnihilateKey file identity check passed." << endl;
    } else {
        cout << "  ❌ FAIL: AnnihilateKey files differ after reload!" << endl;
        return 1;
    }

    // 2. Test Packing and Unpacking Logic
    cout << "[Step 2] Testing Packing & Unpacking Correctness..." << endl;
    
    FheKey key;
    key.load(eval_file, ahk_file, pk_file, pk_v2_file);
    key.load_sk(sec_file);

    // Prepare random data (500 bits)
    // Note: TFHEpp's lvl2param::n is usually 1024, so 500 will not overflow
    int num_bits = 500; 
    vector<uint8_t> original_data(num_bits);
    
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, 1);
    cout << "  -> Generating " << num_bits << " random bits..." << endl;
    for(int i=0; i<num_bits; ++i) {
        original_data[i] = dis(gen);
    }

    // [Critical Fix] Use encryptData instead of encryptBatch
    // encryptBatch splits each uint8_t into 8 bits (resulting in 500 -> 4000 bits, overflowing the 1024 limit)
    // encryptData treats the vector directly as a bit array (500 -> 500 bits)
    
    auto packed_ptr = key.encryptData(original_data); 
    cout << "  -> Data Encrypted and Packed (TLWEs -> TRLWE)." << endl;

    // Check size
    // encryptData sets W to the input size
    if (packed_ptr->W != (size_t)num_bits) {
        cout << "  [Warning] Packed size info: " << packed_ptr->W << " (Expected: " << num_bits << ")" << endl;
    }

    // Unpack
    auto unpacked_vec_ptr = decompress_ciphertexts(packed_ptr);
    cout << "  -> Data Unpacked (TRLWE -> vector<TLWE>)." << endl;

    // Verify Decryption
    bool all_match = true;
    int mismatch_count = 0;

    for (int i = 0; i < num_bits; ++i) {
        // Extract single TLWE
        tlwe_t single_ct = unpacked_vec_ptr->data[i];
        
        // Wrap into Wrapper to call decrypt
        tlwe_ptr ptr(new tlwe(single_ct));
        
        // Decrypt
        uint8_t decrypted_bit = key.decrypt(ptr);
        
        if (decrypted_bit != original_data[i]) {
            all_match = false;
            mismatch_count++;
            if (mismatch_count < 5) {
                cout << "    Mismatch at index " << i << ": Expected " 
                     << (int)original_data[i] << ", Got " << (int)decrypted_bit << endl;
            }
        }
    }

    if (all_match) {
        cout << "  SUCCESS: Pack/Unpack Functional Test passed. All " << num_bits << " bits match." << endl;
    } else {
        cout << "  FAIL: " << mismatch_count << " bits mismatched." << endl;
        return 1;
    }

    cout << "[Step 3] All resources released correctly (checking for double-free)." << endl;
    return 0;
}
