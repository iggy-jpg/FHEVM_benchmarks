#include "fhe_key.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <thread>
#include <cstring>
#include <mutex>
#include <cmath>
using namespace std;

string hex_encode(const vector<uint8_t> &data)
{
    string s = "0x";
    char buf[4];
    for (auto d : data)
    {
        sprintf(buf, "%02x", d);
        s += buf;
    }
    return s;
}

vector<uint8_t> hex_decode(const string &str)
{
    string s = str.substr(2, str.size() - 1);
    vector<uint8_t> bytes;
    bytes.reserve(s.size() / 2);

    for (size_t i = 0; i < s.size(); i += 2)
    {
        unsigned int n;
        sscanf(&s[i], "%2x", &n);
        bytes.push_back(static_cast<uint8_t>(n));
    }
    return bytes;
}

bool is_exist(string dirname, string handle)
{
    std::ifstream file(dirname + "/" + handle, std::ios::binary);
    if (file.is_open())
    {
        file.close();
        return true;
    }
    return false;
}

enum class DataType {
    Bool = 0,
    Uint8 = 1,
    Uint16 = 2,
    Uint32 = 3,
    Uint64 = 4,
    Uint128 = 5,
    Uint256 = 6,
};

string decrypt_ciphertext(FheKey& key, const tlwe_vec_ptr& ciphertext) {
    string result;
    if (ciphertext->W == 1) {
        tlwe_ptr ct = to_tlwe_ptr(ciphertext);
        uint8_t pt = key.decrypt(ct);
        result = pt == 0? "0": "1";
    }
    else {
        auto bytes = key.decryptBatch(ciphertext);
        result = hex_encode(bytes);
    }
    return result;
}

trlwe_ptr trivial_encrypt(FheKey &key, const vector<uint8_t> &bytes, DataType data_type)
{
    int size = 0;
    uint8_t type = (int)data_type;
    if (type == 0)
    {
        size = 1;
    }
    else if (type >= 1 && type <= 6)
    {
        size = 1 << (type + 2);
    }
    else
    {
        throw std::runtime_error("unsupported type");
    }

    trlwe_ptr ct;
    if (size == 1)
    {
        ct = key.trivialEncrypt(bytes[0]);
    }
    else
    {
        size_t num_bytes = size >> 3;
        vector<uint8_t> tmp(bytes.data() + (bytes.size() - num_bytes), bytes.data() + bytes.size());
        ct = key.trivialEncryptBatch(tmp);
    }
    return ct;
}
string plain_add(string lhs, string rhs) {
    vector<uint8_t> lhs_bytes = hex_decode(lhs);
    vector<uint8_t> rhs_bytes = hex_decode(rhs);
    bool carry = false;
    vector<uint8_t> res_bytes;
    res_bytes.resize(lhs_bytes.size());
    for (size_t i = lhs_bytes.size(); i > 0; i--) {
        uint8_t bit = lhs_bytes[i - 1] + rhs_bytes[i - 1];
        res_bytes[i - 1] = carry? bit + 1: bit;
        carry = res_bytes[i - 1] < lhs_bytes[i - 1];
    }
    return hex_encode(res_bytes);
}

string plain_sub(string lhs, string rhs) {
    vector<uint8_t> lhs_bytes = hex_decode(lhs);
    vector<uint8_t> rhs_bytes = hex_decode(rhs);
    bool carry = false;
    vector<uint8_t> res_bytes;
    res_bytes.resize(lhs_bytes.size());
    for (size_t i = lhs_bytes.size(); i > 0; i--) {
        uint8_t bit = lhs_bytes[i - 1] - rhs_bytes[i - 1];
        res_bytes[i - 1] = carry? bit - 1: bit;
        carry = res_bytes[i - 1] > lhs_bytes[i - 1];
    }
    return hex_encode(res_bytes);
}

string plain_lt(string lhs, string rhs) {
    vector<uint8_t> lhs_bytes = hex_decode(lhs);
    vector<uint8_t> rhs_bytes = hex_decode(rhs);
    vector<uint8_t> res_bytes;
    res_bytes.resize(lhs_bytes.size());
    for (size_t i = 0; i < lhs.size(); i++) {
        if (lhs_bytes[i] < rhs_bytes[i]) {
            return "1";
        }
        else if (lhs_bytes[i] > rhs_bytes[i]) {
            return "0";
        }
    }
    return "0";
}

string plain_select(string lhs, string rhs, string ctl) {
    if (ctl == "1") {
        return rhs;
    }
    return lhs;
}

tlwe_vec_ptr trivial_encrypt_integer(FheKey& key, string hex) {
    vector<uint8_t> bytes = hex_decode(hex);
    trlwe_ptr ct = key.encryptBatch(bytes);
    auto ct2 = decompress_ciphertexts(ct);
    return ct2;
}

void run()
{
    FheKey key;
    FheKey decrypt_key;
    map<string, trlwe_ptr> global_ct_map;

    try
    {
        key.load(
            "server_eval.bin",
            "server_ahk.bin",
            "server_pk.bin",
            "server_pk_v2.bin");
        decrypt_key.load_sk("client_secret.bin");
    }
    catch (const std::exception &e)
    {
        std::cerr << "[Thread " << std::this_thread::get_id()
                  << "] ❌ 密钥加载失败: " << e.what() << std::endl;
        return;
    }
    auto ct6 = trivial_encrypt_integer(key, "0x000000000000000000000000000000000000000000000000000000003b9aca01");
    auto ct_upload0 = trivial_encrypt_integer(key, "0x000000000000000000000000000000000000000000000000000000003b9aca00");

    auto ct1 = trivial_encrypt_integer(key, "0x000000000000000000000000000000000000000000000000000000003b9aca00");
    auto ct2 = trivial_encrypt_integer(key, "0x0000000000000000000000000000000000000000000000000000000000000000");
    auto ct3 = key.add_unsigned_pal(ct2, ct1); 
    {
        auto pt1 = decrypt_ciphertext(decrypt_key, ct1);
        auto pt2 = decrypt_ciphertext(decrypt_key, ct2);
        auto pt3 = decrypt_ciphertext(decrypt_key, ct3);
        std::cout << pt2 << " + " << pt1 << " ?= " << pt3 << std::endl;
        auto pt = plain_add(pt2, pt1);
        std::cout << "expected: " << pt << std::endl;
        if (pt != pt3) throw std::runtime_error("compute error");
    }
    auto ct4 = trivial_encrypt_integer(key, "0x0000000000000000000000000000000000000000000000000000000000000000");
    auto ct5 = key.add_unsigned_pal(ct4, ct1); 
    {
        auto pt1 = decrypt_ciphertext(decrypt_key, ct1);
        auto pt4 = decrypt_ciphertext(decrypt_key, ct4);
        auto pt5 = decrypt_ciphertext(decrypt_key, ct5);
        std::cout << pt4 << " + " << pt1 << " ?= " << pt5 << std::endl;
        auto pt = plain_add(pt4, pt1);
        std::cout << "expected: " << pt << std::endl;
        if (pt != pt5) throw std::runtime_error("compute error");
    }
    auto ct7 = key.lt_unsigned_pal(ct5, ct6);
    {
        auto pt5 = decrypt_ciphertext(decrypt_key, ct5);
        auto pt6 = decrypt_ciphertext(decrypt_key, ct6);
        auto pt7 = decrypt_ciphertext(decrypt_key, ct7);
        std::cout << pt5 << " < " << pt6 << " ?= " << pt7 << std::endl;
        auto pt = plain_lt(pt5, pt6);
        std::cout << "expected: " << pt << endl;
        if (pt != pt7) throw std::runtime_error("compute error");
    }
    auto ct8 = key.sub_unsigned_pal(ct5, ct6);
    {
        auto pt5 = decrypt_ciphertext(decrypt_key, ct5);
        auto pt6 = decrypt_ciphertext(decrypt_key, ct6);
        auto pt8 = decrypt_ciphertext(decrypt_key, ct8);
        std::cout << pt5 << " - " << pt6 << " ?= " << pt8 << std::endl;
        auto pt = plain_sub(pt5, pt6);
        std::cout << "expect: " << pt << std::endl;
        if (pt != pt8) throw std::runtime_error("compute error");
    }
    auto ct9 = key.select(ct5, ct8, ct7);
    {
        auto pt5 = decrypt_ciphertext(decrypt_key, ct5);
        auto pt7 = decrypt_ciphertext(decrypt_key, ct7);
        auto pt8 = decrypt_ciphertext(decrypt_key, ct8);
        auto pt9 = decrypt_ciphertext(decrypt_key, ct9);
        std::cout << pt7 << " ? " << pt5 << " : " << pt8 << " ?= " << pt9 << std::endl;
        auto pt = plain_select(pt5, pt8, pt7);
        std::cout << "expected: " << pt << std::endl;
        if (pt != pt9) throw std::runtime_error("compute error");
    }
    auto ct10 = trivial_encrypt_integer(key, "0x0000000000000000000000000000000000000000000000000000000000000000");
    auto ct11 = key.add_unsigned_pal(ct10, ct6);
    {
        auto pt6 = decrypt_ciphertext(decrypt_key, ct6);
        auto pt10 = decrypt_ciphertext(decrypt_key, ct10);
        auto pt11 = decrypt_ciphertext(decrypt_key, ct11);
        std::cout << pt10 << " + " << pt6 << " ?= " << pt11 << std::endl;
        auto pt = plain_add(pt10, pt6);
        std::cout<< "expected: " << pt << std::endl;
        if (pt != pt11) throw std::runtime_error("compute error");
    }
    auto ct12 = trivial_encrypt_integer(key, "0x0000000000000000000000000000000000000000000000000000000000000000");
    auto ct13 = key.select(ct12, ct11, ct7);
    {
        auto pt7 = decrypt_ciphertext(decrypt_key, ct7);
        auto pt11 = decrypt_ciphertext(decrypt_key, ct11);
        auto pt12 = decrypt_ciphertext(decrypt_key, ct12);
        auto pt13 = decrypt_ciphertext(decrypt_key, ct13);
        std::cout << pt7 << " ? " << pt12 << " : " << pt11 << " ?= " << pt13 << std::endl;
        auto pt = plain_select(pt12, pt11, pt7);
        std::cout << "expected: " << pt << endl;
        if (pt != pt13) throw std::runtime_error("compute error");
    }
}

int main(int argc, const char *argv[])
{
    FheKey generator_key;
    std::cout << "[Main] Checking/Initializing Base FHE Keys..." << std::endl;
    
    bool keys_complete = is_exist(".", "client_secret.bin")
                      && is_exist(".", "server_eval.bin")
                      && is_exist(".", "server_ahk.bin")
                      && is_exist(".", "server_pk.bin")
                      && is_exist(".", "server_pk_v2.bin");
    if (!keys_complete)
    {
        std::cout << "[Main] Key files incomplete or not found. Generating..." << std::endl;
        generator_key.keygen(0);
        generator_key.save(
            "client_secret.bin",
            "server_eval.bin",
            "server_ahk.bin",
            "server_pk.bin",
            "server_pk_v2.bin");
    }
    else
    {
        std::cout << "[Main] Base keys exist. Threads will copy them." << std::endl;
    }

    run();

    return 0;
}
