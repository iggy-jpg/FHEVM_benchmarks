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
#include <cstdlib>
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

vector<string> split(string str, string delimiter)
{
    vector<string> result;
    auto beg = 0;
    auto end = str.find(delimiter, beg);
    while (end != string::npos)
    {
        result.push_back(str.substr(beg, end - beg));
        beg = end + 1;
        end = str.find(delimiter, beg);
    }
    if (beg < str.size())
    {
        result.push_back(str.substr(beg, str.size() - beg));
    }
    return result;
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




void load_ciphertexts(map<string, trlwe_ptr> &ct_map, const char *file, FheKey &key)
{
    std::ifstream fp(file, std::ios::in);
    if (!fp.is_open())
    {
        std::cerr << "fail to open file " << file << std::endl;
        return;
    }
    string line;
    while (std::getline(fp, line))
    {
        auto pos = line.find(":");
        string handle = line.substr(0, pos);
        string value = line.substr(pos + 1, line.size() - (pos + 1));
        vector<uint8_t> bytes = hex_decode(value);
        trlwe_ptr ct = key.encryptBatch(bytes);
        ct_map.emplace(handle, std::move(ct));
    }
}

void print_pt(string hint, tlwe_vec_ptr &ct, FheKey &key)
{
     
    // vector<uint8_t> bytes = key.decryptBatch(ct);
    // string pt = hex_encode(bytes);
    // std::lock_guard<std::mutex> lock(io_mutex);
    // std::cout << hint << pt << std::endl;
}

void save_ciphertexts(string handle, tlwe_vec_ptr &ct, FheKey &key)
{
    // print_pt("save ciphertexts ", ct, key);
    std::ofstream file("ciphertexts/" + handle, std::ios::binary);
    auto bytes = serialize_ciphertexts(ct);
    if (bytes.empty())
    {
        throw std::runtime_error("empty bytes");
    }
    
    // print_pt("recover ct ", ct2, key);
    file.write(reinterpret_cast<const char *>(bytes.data()), bytes.size());
    file.close();
}

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

trlwe_ptr trivial_encrypt(FheKey &key, const vector<uint8_t> &bytes, string handle)
{
    int size = 0;
    vector<uint8_t> handle_bytes = hex_decode(handle);
    uint8_t type = handle_bytes[30];
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

vector<tlwe_vec_ptr> get_ciphertexts(string handle, const vector<string> &inputs, int op, int is_scalar, FheKey &key)
{
    vector<tlwe_vec_ptr> result;
    result.reserve(inputs.size());

    for (size_t i = 0; i < inputs.size(); i++)
    {
        if ((i == 1 && is_scalar == 1) || op == 16)
        {
            
            vector<uint8_t> bytes = hex_decode(inputs[i]);
            if (bytes.empty())
            {
                throw std::runtime_error("[Error] hex_decode returned empty bytes for handle: " + handle);
            }

            trlwe_ptr ct;
            if (op == 16)
            {
                ct = trivial_encrypt(key, bytes, handle);
            }
            else
            {
                ct = trivial_encrypt(key, bytes, inputs[0]);
            }

            
            if (!ct)
            {
                throw std::runtime_error("[Error] trivial_encrypt returned nullptr for handle: " + handle);
            }

            auto ct2 = decompress_ciphertexts(ct);
            if (!ct2)
            {
                throw std::runtime_error("[Error] decompress_ciphertexts returned nullptr");
            }

            result.push_back(std::move(ct2));
        }
        else
        {
            
            string filepath = "ciphertexts/" + inputs[i];
            std::ifstream file(filepath, std::ios::binary);

            
            if (!file.is_open())
            {
                throw std::runtime_error("[Error] Missing dependency! Cannot open file: " + filepath);
            }

            file.seekg(0, std::ios::end);
            
            std::streamsize size = file.tellg();
            file.seekg(0, std::ios::beg);

            
            if (size <= 0)
            {
                throw std::runtime_error("[Error] File is empty or corrupted: " + filepath + " (size: " + to_string(size) + ")");
            }

            vector<char> data(size);

            
            if (!file.read(data.data(), size))
            {
                throw std::runtime_error("[Error] Failed to read complete data from file: " + filepath);
            }

            
            vector<uint8_t> bytes(data.begin(), data.end());

            auto ct = deserialize_ciphertexts(bytes);

            
            if (!ct)
            {
                throw std::runtime_error("[Error] deserialize_ciphertexts returned nullptr for file: " + filepath);
            }

            result.push_back(std::move(ct));
        }
    }
    return result;
}





static std::map<string, Bytes> g_plain_values;
static std::mutex g_plain_mu;

static void record_plain_value(const string& handle, const Bytes& bytes)
{
    std::lock_guard<std::mutex> lk(g_plain_mu);
    g_plain_values[handle] = bytes;
}

static bool lookup_plain_value(const string& handle, Bytes& out)
{
    std::lock_guard<std::mutex> lk(g_plain_mu);
    auto it = g_plain_values.find(handle);
    if (it == g_plain_values.end()) return false;
    out = it->second;
    return true;
}

tlwe_vec_ptr compute(int op, const vector<tlwe_vec_ptr> &cts,
                     const vector<std::unique_ptr<Bytes>> &plains, FheKey &key)
{
    tlwe_vec_ptr res;
    if (op == 1)
    {
        res = key.add_unsigned_pal(cts[0], cts[1]);
    }
    else if (op == 2)
    {
        res = key.sub_unsigned_pal(cts[0], cts[1]);
    }
    else if (op == 5)
    {
        if (plains.size() > 1 && plains[1])      res = key.and_plain(cts[0], *plains[1]);
        else if (plains.size() > 0 && plains[0]) res = key.and_plain(cts[1], *plains[0]);
        else res = key.bit_and(cts[0], cts[1]);
    }
    else if (op == 6)
    {
        if (plains.size() > 1 && plains[1])      res = key.or_plain(cts[0], *plains[1]);
        else if (plains.size() > 0 && plains[0]) res = key.or_plain(cts[1], *plains[0]);
        else res = key.bit_or(cts[0], cts[1]);
    }
    else if (op == 8)
    {
        res = key.bit_not(cts[0]);
    }
    else if (op == 11)
    {
        if (plains.size() > 1 && plains[1])      res = key.ge_plain_pal(cts[0], *plains[1]);
        else if (plains.size() > 0 && plains[0]) res = key.le_plain_pal(cts[1], *plains[0]);
        else res = key.ge_unsigned_pal(cts[0], cts[1]);
    }
    else if (op == 14)
    {
        if (plains.size() > 1 && plains[1])      res = key.lt_plain_pal(cts[0], *plains[1]);
        else if (plains.size() > 0 && plains[0]) res = key.gt_plain_pal(cts[1], *plains[0]);
        else res = key.lt_unsigned_pal(cts[0], cts[1]);
    }
    else if (op == 15)
    {
        
        
        if (plains.size() > 2 && plains[1] && plains[2] && cts[0]->W >= 1) {
            
            res = key.select_plain(*plains[2], *plains[1], cts[0], cts[1]->W);
        } else {
            res = key.select(cts[2], cts[1], cts[0]);
        }
    }
    else if (op == 16)
    {
        res.reset(new tlwe_vec(*cts[0].get()));
    }
    else
    {
        throw std::runtime_error("not supported op");
    }
    auto pt = decrypt_ciphertext(key, res);
     







    auto pt2 = decrypt_ciphertext(key, res);
    return res;
}



tlwe_vec_ptr run_fhe_op(string handle, int op, string deps, int is_scalar, const map<string, trlwe_ptr> &ct_map, FheKey &key)
{
    vector<string> inputs = split(deps, ";");
    if (op == 17)
    {
        auto iter = ct_map.find(handle);
        
        auto ct = decompress_ciphertexts(iter->second);
        save_ciphertexts(handle, ct, key);
        return ct;
    }

    
    if (op == 16 && !inputs.empty())
    {
        record_plain_value(handle, hex_decode(inputs[0]));
    }
    
    static const bool plain_fold_enabled =
        (std::getenv("AUTOHOG_DISABLE_PLAIN_FOLD") == nullptr);
    vector<std::unique_ptr<Bytes>> plains(inputs.size());
    for (size_t i = 0; plain_fold_enabled && i < inputs.size(); i++)
    {
        if (i == 1 && is_scalar == 1)
        {
            plains[i].reset(new Bytes(hex_decode(inputs[i])));
        }
        else
        {
            Bytes v;
            if (lookup_plain_value(inputs[i], v)) plains[i].reset(new Bytes(v));
        }
    }

    
    vector<tlwe_vec_ptr> cts = get_ciphertexts(handle, inputs, op, is_scalar, key);

    
    auto res = compute(op, cts, plains, key);

    save_ciphertexts(handle, res, key);
    return res;
}

void run()
{
    FheKey local_key;
    FheKey decrypt_key;
    map<string, trlwe_ptr> global_ct_map;

    try
    {
        
        
        local_key.load(
            "server_eval.bin",
            "server_ahk.bin",
            "server_pk.bin",
            "server_pk_v2.bin",
            "zk_crs.bin");
        decrypt_key.load_sk("client_secret.bin");
    }
    catch (const std::exception &e)
    {
        std::cerr << "[run] ❌ 密钥加载失败: " << e.what() << std::endl;
        return;
    }

    // [Sanity Check]
    {
        tlwe_vec_t bits1(32, TFHEpp::CT_ZERO); bits1[0] = TFHEpp::CT_ONE;
        tlwe_vec_ptr ct1(new tlwe_vec(bits1));

        tlwe_vec_t bits0(32, TFHEpp::CT_ZERO);
        tlwe_vec_ptr ct0(new tlwe_vec(bits0));

        tlwe_vec_t s1_bits = {TFHEpp::CT_ONE};
        tlwe_vec_ptr sel1(new tlwe_vec(s1_bits));

        tlwe_vec_t s0_bits = {TFHEpp::CT_ZERO};
        tlwe_vec_ptr sel0(new tlwe_vec(s0_bits));

        // FheKey::select(lhs, rhs, sel) = sel ? rhs : lhs；
        
        auto res110 = local_key.select(ct1, ct0, sel1);
        auto pt110 = decrypt_key.decryptBatch(res110);
        cout << "[Sanity Check] select(lhs=1, rhs=0, sel=1) = " << (int)pt110[0] << " (expected 0 = rhs)" << endl;

        auto res010 = local_key.select(ct1, ct0, sel0);
        auto pt010 = decrypt_key.decryptBatch(res010);
        cout << "[Sanity Check] select(lhs=1, rhs=0, sel=0) = " << (int)pt010[0] << " (expected 128 = lhs MSB)" << endl;
    }



    load_ciphertexts(global_ct_map, "data/ciphertexts_upload_select.txt", local_key);

    map<string, trlwe_ptr> ct_map;

    {
        std::ifstream fp("data/computations_select.txt", std::ios::in);
        string line;
        while (std::getline(fp, line))
        {
            vector<string> strs = split(line, ":");
            string handle = strs[0];
            int op = atoi(strs[1].c_str());
            string deps = strs[2];
            int is_scalar = atoi(strs[3].c_str());

            
            auto ciphertext = run_fhe_op(handle, op, deps, is_scalar, global_ct_map, local_key);
            auto plaintext = decrypt_ciphertext(decrypt_key, ciphertext);

            {
                cout << "execute handle: " << handle << "  op: " << std::setw(2) << op <<  "  plaintext: " << plaintext << endl;
            }
        }

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
                      && is_exist(".", "server_pk_v2.bin")
                      && is_exist(".", "zk_crs.bin");
    if (!keys_complete)
    {
        std::cout << "[Main] Key files incomplete or not found. Generating..." << std::endl;
        generator_key.keygen(0);
        generator_key.save(
            "client_secret.bin",
            "server_eval.bin",
            "server_ahk.bin",
            "server_pk.bin",
            "server_pk_v2.bin",
            "zk_crs.bin");
    }
    else
    {
        std::cout << "[Main] Base keys exist. Threads will copy them." << std::endl;
    }

    run();

    return 0;
}
