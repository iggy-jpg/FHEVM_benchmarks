#include "fhe_key.h"
#include <algorithm>
#include <sys/stat.h>
#include <sys/types.h>
#if defined(_WIN32)
#include <direct.h>
#endif
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
#include <atomic>
#include <cerrno>
#include <cfenv>
#include <chrono>
#include <filesystem>
#include <iomanip>
#include <iterator>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <openssl/evp.h>
#if defined(__SSE2__)
#include <immintrin.h>
#endif
using namespace std;

std::mutex io_mutex;
std::atomic<bool> worker_failed{false};
std::mutex worker_error_mutex;
std::vector<std::string> worker_errors;
std::string expected_constants_hash;
std::string consensus_work_root;
int operator_threads_per_worker = 1;

static std::string worker_directory(int index)
{
    return consensus_work_root + "/worker_" + std::to_string(index);
}

class ConsensusWorkDirectory
{
  public:
    explicit ConsensusWorkDirectory(int workers) : workers_(workers)
    {
        const char *configured = std::getenv("AUTOHOG_CONSENSUS_WORKDIR");
        if (configured != nullptr && *configured != '\0') {
            consensus_work_root = configured;
        } else {
            const auto nonce = std::chrono::high_resolution_clock::now()
                                   .time_since_epoch().count();
            consensus_work_root =
                (std::filesystem::temp_directory_path() /
                 ("autohog_consensus_" + std::to_string(nonce))).string();
        }

        std::error_code error;
        root_created_ = std::filesystem::create_directories(consensus_work_root, error);
        if (error)
            throw std::runtime_error("cannot create consensus work directory: " +
                                     consensus_work_root + ": " + error.message());

        // A caller-provided work root may be reused. Remove only this test's
        // precisely named worker directories, never unrelated files.
        for (int i = 1; i <= workers_; ++i) {
            std::filesystem::remove_all(worker_directory(i), error);
            if (error)
                throw std::runtime_error("cannot reset " + worker_directory(i) +
                                         ": " + error.message());
        }
        std::cout << "[Main] work_directory=" << consensus_work_root << std::endl;
    }

    ~ConsensusWorkDirectory()
    {
        if (std::getenv("AUTOHOG_KEEP_CONSENSUS_WORKDIR") != nullptr) return;
        std::error_code error;
        for (int i = 1; i <= workers_; ++i) {
            std::filesystem::remove_all(worker_directory(i), error);
            if (error) {
                std::cerr << "[Cleanup] cannot remove " << worker_directory(i)
                          << ": " << error.message() << std::endl;
                error.clear();
            }
        }
        if (root_created_) std::filesystem::remove(consensus_work_root, error);
    }

  private:
    int workers_;
    bool root_created_ = false;
};

static const std::string &computations_file()
{
    static const std::string path = [] {
        const char *env = std::getenv("AUTOHOG_CONSENSUS_COMPUTATIONS");
        return std::string(env ? env : "data/computations.txt");
    }();
    return path;
}

static int make_directory(const std::string &path)
{
#if defined(_WIN32)
    return _mkdir(path.c_str());
#else
    return mkdir(path.c_str(), 0755);
#endif
}

static std::string sha256_bytes(const uint8_t *data, size_t size)
{
    EVP_MD_CTX *ctx = EVP_MD_CTX_new();
    if (ctx == nullptr) throw std::runtime_error("EVP_MD_CTX_new failed");
    unsigned char digest[EVP_MAX_MD_SIZE];
    unsigned int digest_size = 0;
    const bool ok = EVP_DigestInit_ex(ctx, EVP_sha256(), nullptr) == 1 &&
                    EVP_DigestUpdate(ctx, data, size) == 1 &&
                    EVP_DigestFinal_ex(ctx, digest, &digest_size) == 1;
    EVP_MD_CTX_free(ctx);
    if (!ok) throw std::runtime_error("SHA-256 calculation failed");

    std::ostringstream out;
    out << std::hex << std::setfill('0');
    for (unsigned int i = 0; i < digest_size; ++i)
        out << std::setw(2) << static_cast<unsigned int>(digest[i]);
    return out.str();
}

static std::string sha256_bytes(const Bytes &bytes)
{
    return sha256_bytes(bytes.data(), bytes.size());
}

static std::string constants_sha256()
{
    Bytes bytes(sizeof(CT_ZERO) + sizeof(CT_ONE));
    std::memcpy(bytes.data(), CT_ZERO.data(), sizeof(CT_ZERO));
    std::memcpy(bytes.data() + sizeof(CT_ZERO), CT_ONE.data(), sizeof(CT_ONE));
    return sha256_bytes(bytes);
}

static unsigned int configure_floating_point_environment()
{
    if (std::fesetround(FE_TONEAREST) != 0)
        throw std::runtime_error("failed to set FE_TONEAREST");
#if defined(__SSE2__)
    unsigned int csr = _mm_getcsr();
    csr &= ~static_cast<unsigned int>(_MM_ROUND_MASK);
    csr |= _MM_ROUND_NEAREST;
    csr |= _MM_FLUSH_ZERO_ON;
    csr |= (1U << 6); // DAZ
    _mm_setcsr(csr);
    return _mm_getcsr();
#else
    return 0;
#endif
}

// ==========================================

// ==========================================


bool copy_file(const std::string &src, const std::string &dst)
{
    std::ifstream source(src, std::ios::binary);
    std::ofstream dest(dst, std::ios::binary);
    if (!source || !dest)
        return false;
    dest << source.rdbuf();
    return true;
}

static std::string calculate_file_hash(const std::string &filename)
{
    return ZKProof::compute_file_sha256_hex(filename);
}



bool compare(const vector<uint8_t> &lhs, const vector<uint8_t> &rhs)
{
    bool b = lhs.size() == rhs.size() && memcmp(lhs.data(), rhs.data(), lhs.size()) == 0;
    return b;
}

string hex_encode(const vector<uint8_t> &data)
{
    string s;
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
        throw std::runtime_error("fail to open file " + std::string(file));

    const char *cache_env = std::getenv("AUTOHOG_CONSENSUS_INPUTS");
    const string cache_dir = cache_env ? cache_env : "";
    if (!cache_dir.empty() && make_directory(cache_dir) != 0 && errno != EEXIST)
        throw std::runtime_error("cannot create input cache " + cache_dir);

    string line;
    while (std::getline(fp, line))
    {
        auto pos = line.find(":");
        if (pos == string::npos)
            throw std::runtime_error("invalid ciphertext upload line");
        string handle = line.substr(0, pos);
        string value = line.substr(pos + 1, line.size() - (pos + 1));
        trlwe_ptr ct;
        const string cache_path = cache_dir.empty() ? "" : cache_dir + "/" + handle;
        if (!cache_path.empty() && is_exist(cache_dir, handle)) {
            ct = load_compressed_ciphertexts(cache_path);
            if (!ct) throw std::runtime_error("cannot load cached input " + cache_path);
        } else {
            vector<uint8_t> bytes = hex_decode(value);
            ct = key.encryptBatch(bytes);
            if (!ct) throw std::runtime_error("cannot encrypt input " + handle);
            if (!cache_path.empty()) save_compressed_ciphertexts(cache_path, ct);
        }
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

Bytes serialize_sanitized_copy(const tlwe_vec_ptr &ct);

void save_ciphertexts(string dirname, string handle, tlwe_vec_ptr &ct, FheKey &key)
{
    // print_pt("save ciphertexts ", ct, key);
    std::ofstream file(dirname + "/" + handle, std::ios::binary);
    auto bytes = serialize_ciphertexts(ct);
    if (bytes.empty())
    {
        throw std::runtime_error("empty bytes");
    }
    
    // print_pt("recover ct ", ct2, key);
    file.write(reinterpret_cast<const char *>(bytes.data()), bytes.size());
    file.close();

    const string diagnostic_dir = dirname + "/sanitized";
    make_directory(diagnostic_dir);
    const Bytes sanitized = serialize_sanitized_copy(ct);
    std::ofstream diagnostic(diagnostic_dir + "/" + handle, std::ios::binary);
    if (!diagnostic)
        throw std::runtime_error("cannot create sanitized diagnostic for " + handle);
    diagnostic.write(reinterpret_cast<const char *>(sanitized.data()), sanitized.size());
    if (!diagnostic)
        throw std::runtime_error("cannot write sanitized diagnostic for " + handle);
}

Bytes serialize_sanitized_copy(const tlwe_vec_ptr &ct)
{
    tlwe_vec_ptr copy(new tlwe_vec(*ct.get()));
    for (auto &elem : copy->data)
    {
        sanitize_tlwe(elem);
    }
    return serialize_ciphertexts(copy);
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

vector<tlwe_vec_ptr> get_ciphertexts(string handle, string dirname, const vector<string> &inputs, int op, int is_scalar, FheKey &key)
{
    static std::mutex fixed_constant_cache_mutex;
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

            std::unique_lock<std::mutex> cache_lock;
            string constant_cache_path;
            if (op == 16) {
                const char *cache_env =
                    std::getenv("AUTOHOG_CONSENSUS_INPUTS");
                if (cache_env != nullptr && *cache_env != '\0') {
                    constant_cache_path = std::string(cache_env) +
                                          "/constant_" + handle;
                    cache_lock = std::unique_lock<std::mutex>(
                        fixed_constant_cache_mutex);
                    std::ifstream cached(constant_cache_path,
                                         std::ios::binary);
                    if (cached) {
                        Bytes serialized((std::istreambuf_iterator<char>(cached)),
                                         std::istreambuf_iterator<char>());
                        auto cached_ct = deserialize_ciphertexts(serialized);
                        if (!cached_ct)
                            throw std::runtime_error(
                                "cannot deserialize cached constant " +
                                constant_cache_path);
                        result.push_back(std::move(cached_ct));
                        continue;
                    }
                }
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

            if (!constant_cache_path.empty()) {
                const Bytes serialized = serialize_ciphertexts(ct2);
                std::ofstream cached(constant_cache_path, std::ios::binary);
                if (!cached || !cached.write(
                        reinterpret_cast<const char *>(serialized.data()),
                        serialized.size()))
                    throw std::runtime_error(
                        "cannot write cached constant " + constant_cache_path);
            }
            result.push_back(std::move(ct2));
        }
        else
        {
            
            string filepath = dirname + "/" + inputs[i];
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
        res = key.add(cts[0], cts[1]);
    }
    else if (op == 2)
    {
        res = key.sub(cts[0], cts[1]);
    }
    else if (op == 3)
    {
        // Keep the parallel circuit topology identical at every external
        // concurrency; only its scheduling width changes.
        res = key.mul_pal(cts[0], cts[1], operator_threads_per_worker);
    }
    else if (op == 4)
    {
        res = key.div_pal(cts[0], cts[1], operator_threads_per_worker);
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
        
        if (plains.size() > 1 && plains[1])      res = key.ge_plain_signed(cts[0], *plains[1]);
        else if (plains.size() > 0 && plains[0]) res = key.le_plain_signed(cts[1], *plains[0]);
        else res = key.ge(cts[0], cts[1]);
    }
    else if (op == 14)
    {
        if (plains.size() > 1 && plains[1])      res = key.lt_plain_signed(cts[0], *plains[1]);
        else if (plains.size() > 0 && plains[0]) res = key.gt_plain_signed(cts[1], *plains[0]);
        else res = key.lt(cts[0], cts[1]);
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
    if (!res)
        throw std::runtime_error("operator returned a null ciphertext");
    return res;
}



void run_fhe_op(vector<string> &error_handles, string dirname, string handle, int op, string deps, int is_scalar, const map<string, trlwe_ptr> &ct_map, FheKey &key)
{
    (void)error_handles;
    if (is_exist(dirname, handle))
    {
        return;
    }
    vector<string> inputs = split(deps, ";");
    if (op == 17)
    {
        auto iter = ct_map.find(handle);
        if (iter == ct_map.end())
            throw std::runtime_error("missing uploaded ciphertext for " + handle);
        
        auto ct = decompress_ciphertexts(iter->second);
        save_ciphertexts(dirname, handle, ct, key);
        return;
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

    
    vector<tlwe_vec_ptr> cts = get_ciphertexts(handle, dirname, inputs, op, is_scalar, key);

    // Lightweight consensus mode: compute once per thread and compare results
    // only after all threads finish.
    auto res = compute(op, cts, plains, key);
    save_ciphertexts(dirname, handle, res, key);
}

void run(int index, const map<string, trlwe_ptr> &global_ct_map, FheKey &shared_key)
{
    const unsigned int mxcsr = configure_floating_point_environment();
    string dirname = worker_directory(index);
    if (make_directory(dirname) != 0 && errno != EEXIST)
        throw std::runtime_error("cannot create " + dirname);

    
    static const bool use_shared_key =
        std::getenv("AUTOHOG_CONSENSUS_SHARED_KEY") != nullptr;
    vector<string> key_files = {
        "client_secret.bin",
        "server_eval.bin",
        "server_annihilate.bin", 
        "server_pk.bin",
        "server_pk_v2.bin"};

    
    for (const auto &kf : key_files)
    {
        if (use_shared_key) break;
        string dest_path = dirname + "/" + kf;

        
        if (!copy_file(kf, dest_path))
        {
            throw std::runtime_error("failed to copy key file " + kf);
        }

        
        const string src_hash = calculate_file_hash(kf);
        const string dst_hash = calculate_file_hash(dest_path);

        if (src_hash.empty() || src_hash != dst_hash)
        {
            throw std::runtime_error("SHA-256 mismatch after copying " + kf);
        }
    }

    
    FheKey local_key;
    FheKey *worker_key = &shared_key;
    if (!use_shared_key) {
        local_key.load(
            dirname + "/server_eval.bin",
            dirname + "/server_annihilate.bin",
            dirname + "/server_pk.bin",
            dirname + "/server_pk_v2.bin");
        worker_key = &local_key;
    }
    if (constants_sha256() != expected_constants_hash)
        throw std::runtime_error("CT_ZERO/CT_ONE changed while loading worker key");

    {
        std::lock_guard<std::mutex> lock(io_mutex);
        std::cout << "[Thread " << index << "] mxcsr=0x" << std::hex << mxcsr
                  << std::dec << " constants_sha256=" << expected_constants_hash << std::endl;
    }

    map<string, trlwe_ptr> ct_map;
    
    // load_ciphertexts(ct_map, "data/ciphertexts_upload.txt", local_key);

    vector<string> error_handles;
    {
        std::ifstream fp(computations_file(), std::ios::in);
        if (!fp) throw std::runtime_error("cannot open " + computations_file());
        string line;
        while (std::getline(fp, line))
        {
            if (line.empty() || line[0] == '#') continue;
            vector<string> strs = split(line, ":");
            if (strs.size() < 4)
                throw std::runtime_error("invalid computations.txt line");
            string handle = strs[0];
            int op = atoi(strs[1].c_str());
            string deps = strs[2];
            int is_scalar = atoi(strs[3].c_str());

            
            run_fhe_op(error_handles, dirname, handle, op, deps, is_scalar,
                       global_ct_map, *worker_key);

        }

        if (!error_handles.empty())
        {
            std::lock_guard<std::mutex> lock(io_mutex);
            cout << "=============" << dirname << "===========" << endl;
            for (auto &h : error_handles)
            {
                cout << h << endl;
            }
        }
    }
}

// void verify_consistency(int num_threads) {
//     if (num_threads < 2) {

//         return;
//     }

//     std::cout << "\n=========================================" << std::endl;


//     std::ifstream fp("data/computations.txt", std::ios::in);
//     if (!fp.is_open()) {

//         return;
//     }

//     std::string line;
//     int match_count = 0;
//     int mismatch_count = 0;

//     while (std::getline(fp, line)) {
//         std::vector<std::string> strs = split(line, ":");
//         if (strs.empty()) continue;
//         std::string handle = strs[0];


//         std::string base_file = "tmp_1/" + handle;
//         std::ifstream bf(base_file, std::ios::binary | std::ios::ate);


//         std::streamsize size = bf.tellg();
//         bf.seekg(0, std::ios::beg);
//         std::vector<char> base_data(size);
//         if (!bf.read(base_data.data(), size)) continue;

//         bool is_consistent = true;

//         for (int i = 2; i <= num_threads; ++i) {
//             std::string cmp_file = "tmp_" + std::to_string(i) + "/" + handle;
//             std::ifstream cf(cmp_file, std::ios::binary | std::ios::ate);

//             if (!cf.is_open()) {
//                 is_consistent = false;

//                 break;
//             }

//             std::streamsize csize = cf.tellg();
//             if (csize != size) {
//                 is_consistent = false;
//                 break;
//             }

//             cf.seekg(0, std::ios::beg);
//             std::vector<char> cmp_data(csize);
//             if (!cf.read(cmp_data.data(), csize)) continue;


//             if (std::memcmp(base_data.data(), cmp_data.data(), size) != 0) {
//                 is_consistent = false;
//                 break;
//             }
//         }

//         if (is_consistent) {
//             match_count++;
//         } else {
//             mismatch_count++;

//         }
//     }




//     if (mismatch_count == 0 && match_count > 0) {

//     } else {

//     }
//     std::cout << "=========================================\n" << std::endl;
// }

static Bytes read_binary_file(const string &path)
{
    std::ifstream file(path, std::ios::binary | std::ios::ate);
    if (!file) throw std::runtime_error("missing result file: " + path);
    const std::streamsize size = file.tellg();
    if (size <= 0) throw std::runtime_error("empty result file: " + path);
    file.seekg(0, std::ios::beg);
    Bytes bytes(static_cast<size_t>(size));
    if (!file.read(reinterpret_cast<char *>(bytes.data()), size))
        throw std::runtime_error("failed to read result file: " + path);
    return bytes;
}

struct CoefficientDiff {
    bool shape_match = true;
    size_t count = 0;
    size_t above_low4 = 0;
    size_t first = std::numeric_limits<size_t>::max();
    int highest_bit = -1;
};

static CoefficientDiff compare_coefficients(const tlwe_vec_ptr &lhs,
                                            const tlwe_vec_ptr &rhs)
{
    CoefficientDiff diff;
    if (!lhs || !rhs || lhs->W != rhs->W || lhs->data.size() != rhs->data.size()) {
        diff.shape_match = false;
        return diff;
    }
    size_t flat = 0;
    for (size_t i = 0; i < lhs->data.size(); ++i) {
        if (lhs->data[i].size() != rhs->data[i].size()) {
            diff.shape_match = false;
            return diff;
        }
        for (size_t j = 0; j < lhs->data[i].size(); ++j, ++flat) {
            const uint64_t delta = static_cast<uint64_t>(lhs->data[i][j]) ^
                                   static_cast<uint64_t>(rhs->data[i][j]);
            if (delta == 0) continue;
            if (diff.count == 0) diff.first = flat;
            ++diff.count;
            if ((delta & ~uint64_t{0xF}) != 0) ++diff.above_low4;
            const int bit = 63 - __builtin_clzll(delta);
            diff.highest_bit = std::max(diff.highest_bit, bit);
        }
    }
    return diff;
}

static map<string, Bytes> load_expected_plaintexts()
{
    map<string, Bytes> expected;
    std::ifstream file("data/consensus_expected.txt");
    if (!file) throw std::runtime_error("cannot open data/consensus_expected.txt");
    string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split(line, ":");
        if (fields.size() != 2)
            throw std::runtime_error("invalid consensus_expected.txt line");
        expected[fields[0]] = hex_decode(fields[1]);
    }
    return expected;
}

static bool compare_result_pair(const string &scope, int num_threads, int op,
                                const string &handle, const string &peer,
                                const string &lhs_path, const string &rhs_path,
                                FheKey &decrypt_key, std::ofstream &summary)
{
    const Bytes lhs_raw = read_binary_file(lhs_path);
    const Bytes rhs_raw = read_binary_file(rhs_path);
    const auto lhs = deserialize_ciphertexts(lhs_raw);
    const auto rhs = deserialize_ciphertexts(rhs_raw);
    if (!lhs || !rhs) throw std::runtime_error("failed to deserialize " + handle);

    const CoefficientDiff diff = compare_coefficients(lhs, rhs);
    const bool raw_match = diff.shape_match && diff.count == 0 && lhs_raw == rhs_raw;
    const Bytes lhs_sanitized = serialize_sanitized_copy(lhs);
    const Bytes rhs_sanitized = serialize_sanitized_copy(rhs);
    const bool sanitized_match = lhs_sanitized == rhs_sanitized;
    const bool plaintext_match = decrypt_key.decryptBatch(lhs) == decrypt_key.decryptBatch(rhs);

    summary << scope << ',' << num_threads << ',' << op << ',' << handle << ',' << peer
            << ',' << (raw_match ? 1 : 0)
            << ',' << (sanitized_match ? 1 : 0)
            << ',' << (plaintext_match ? 1 : 0)
            << ',' << (diff.shape_match ? diff.count : size_t(-1))
            << ',' << diff.above_low4
            << ',' << (diff.first == std::numeric_limits<size_t>::max() ? -1LL
                                                                          : static_cast<long long>(diff.first))
            << ',' << diff.highest_bit
            << ',' << sha256_bytes(lhs_raw) << ',' << sha256_bytes(rhs_raw) << '\n';
    summary.flush();

    if (!raw_match) {
        std::cerr << "[RAW_MISMATCH] scope=" << scope << " threads=" << num_threads
                  << " op=" << op << " handle=" << handle << " peer=" << peer
                  << " coeff_diffs=" << diff.count << " above_low4=" << diff.above_low4
                  << " first=" << diff.first << " highest_bit=" << diff.highest_bit
                  << " sanitized_match=" << sanitized_match
                  << " plaintext_match=" << plaintext_match << std::endl;
    }
    return raw_match && plaintext_match;
}

bool verify_consistency(int num_threads, FheKey &decrypt_key)
{
    const char *summary_env = std::getenv("AUTOHOG_CONSENSUS_SUMMARY");
    const string summary_path = summary_env ? summary_env : "consensus_summary.csv";
    std::ofstream summary(summary_path);
    if (!summary) throw std::runtime_error("cannot create summary: " + summary_path);
    summary << "scope,threads,op,handle,peer,raw_match,sanitized_match,plaintext_match,"
               "coefficient_differences,differences_above_low4,first_difference,"
               "highest_difference_bit,lhs_sha256,rhs_sha256\n";

    const auto expected = load_expected_plaintexts();
    const char *baseline_env = std::getenv("AUTOHOG_CONSENSUS_BASELINE");
    const string baseline_dir = baseline_env ? baseline_env : "";
    if (num_threads == 1 && !baseline_dir.empty()) {
        make_directory(baseline_dir);
        make_directory(baseline_dir + "/sanitized");
    }

    bool all_ok = true;
    int checked = 0;
    std::ifstream computations(computations_file());
    if (!computations) throw std::runtime_error("cannot open " + computations_file());
    string line;
    while (std::getline(computations, line)) {
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split(line, ":");
        if (fields.size() < 4) throw std::runtime_error("invalid computations.txt line");
        const string &handle = fields[0];
        const int op = std::atoi(fields[1].c_str());
        const string current = worker_directory(1) + "/" + handle;
        const Bytes current_raw = read_binary_file(current);
        const auto current_ct = deserialize_ciphertexts(current_raw);
        if (!current_ct) throw std::runtime_error("cannot deserialize " + current);

        const auto expected_it = expected.find(handle);
        if (expected_it != expected.end()) {
            const Bytes actual = decrypt_key.decryptBatch(current_ct);
            const bool expected_ok = actual == expected_it->second;
            std::cout << "[EXPECTED] op=" << op << " handle=" << handle
                      << " result=" << (expected_ok ? "ok" : "incorrect");
            if (!expected_ok) {
                std::cout << " actual=" << hex_encode(actual)
                          << " expected=" << hex_encode(expected_it->second);
            }
            std::cout << std::endl;
            all_ok = expected_ok && all_ok;
        }

        for (int i = 2; i <= num_threads; ++i) {
            all_ok = compare_result_pair("within", num_threads, op, handle,
                                         "worker_" + std::to_string(i), current,
                                         worker_directory(i) + "/" + handle,
                                         decrypt_key, summary) && all_ok;
            ++checked;
        }

        if (!baseline_dir.empty()) {
            const string baseline = baseline_dir + "/" + handle;
            if (num_threads == 1) {
                if (!copy_file(current, baseline))
                    throw std::runtime_error("cannot create baseline for " + handle);
                if (!copy_file(worker_directory(1) + "/sanitized/" + handle,
                               baseline_dir + "/sanitized/" + handle))
                    throw std::runtime_error("cannot create sanitized baseline for " + handle);
            } else {
                all_ok = compare_result_pair("baseline", num_threads, op, handle,
                                             "concurrency_1", baseline, current,
                                             decrypt_key, summary) && all_ok;
                ++checked;
            }
        }
    }

    std::cout << "[Verify] comparisons=" << checked
              << " raw_consensus=" << (all_ok ? "PASS" : "FAIL") << std::endl;
    return all_ok;
}


int main(int argc, const char *argv[])
{
    try {
        int num_threads = 3;
        if (argc > 1) {
            num_threads = std::atoi(argv[1]);
            if (num_threads <= 0) num_threads = 1;
        }
        const char *budget_env =
            std::getenv("AUTOHOG_CONSENSUS_TOTAL_THREADS");
        const unsigned detected = std::thread::hardware_concurrency();
        int total_thread_budget = budget_env != nullptr
            ? std::atoi(budget_env)
            : std::min(16, std::max(1, static_cast<int>(detected)));
        if (total_thread_budget <= 0) total_thread_budget = 1;
        operator_threads_per_worker =
            std::max(1, total_thread_budget / num_threads);
        std::cout << "[Main] external_workers=" << num_threads
                  << " operator_threads_per_worker="
                  << operator_threads_per_worker
                  << " total_thread_budget=" << total_thread_budget
                  << std::endl;
        ConsensusWorkDirectory work_directory(num_threads);

        const unsigned int main_mxcsr = configure_floating_point_environment();
        FheKey generator_key;
        std::cout << "[Main] Checking/Initializing Base FHE Keys..." << std::endl;
        const vector<string> key_files = {
            "client_secret.bin", "server_eval.bin", "server_annihilate.bin",
            "server_pk.bin", "server_pk_v2.bin"};
        bool have_all_keys = true;
        for (const auto &key_file : key_files)
            have_all_keys = have_all_keys && is_exist(".", key_file);

        if (!have_all_keys) {
            std::cout << "[Main] Complete key set not found. Generating..." << std::endl;
            generator_key.keygen(0);
            generator_key.save(
                "client_secret.bin", "server_eval.bin", "server_annihilate.bin",
                "server_pk.bin", "server_pk_v2.bin");
        } else {
            std::cout << "[Main] Loading complete existing key set..." << std::endl;
            generator_key.load(
                "server_eval.bin", "server_annihilate.bin",
                "server_pk.bin", "server_pk_v2.bin");
            generator_key.load_sk("client_secret.bin");
        }

        expected_constants_hash = constants_sha256();
        std::cout << "[Main] mxcsr=0x" << std::hex << main_mxcsr << std::dec
                  << " constants_sha256=" << expected_constants_hash << std::endl;

        map<string, trlwe_ptr> global_ct_map;
        load_ciphertexts(global_ct_map, "data/ciphertexts_upload.txt", generator_key);

        worker_failed.store(false);
        worker_errors.clear();
        std::vector<std::thread> threads;
        for (int i = 1; i <= num_threads; ++i) {
            threads.emplace_back([i, &global_ct_map, &generator_key]() {
                try {
                    run(i, global_ct_map, generator_key);
                } catch (const std::exception &e) {
                    worker_failed.store(true);
                    std::lock_guard<std::mutex> lock(worker_error_mutex);
                    worker_errors.push_back("thread " + std::to_string(i) + ": " + e.what());
                } catch (...) {
                    worker_failed.store(true);
                    std::lock_guard<std::mutex> lock(worker_error_mutex);
                    worker_errors.push_back("thread " + std::to_string(i) + ": unknown exception");
                }
            });
        }
        for (auto &thread : threads) thread.join();

        if (worker_failed.load()) {
            for (const auto &error : worker_errors) std::cerr << "[WorkerError] " << error << '\n';
            return 2;
        }

        std::cout << "[Main] All " << num_threads << " threads done computing." << std::endl;
        return verify_consistency(num_threads, generator_key) ? 0 : 3;
    } catch (const std::exception &e) {
        std::cerr << "[Fatal] " << e.what() << std::endl;
        return 1;
    }
}

// int main(int argc, const char *argv[])
// {

//     int process_id = 1;
//     if (argc > 1)
//     {
//         process_id = std::atoi(argv[1]);
//     }

//     FheKey generator_key;

//     // ==========================================

//     // ==========================================
//     if (process_id == 0)
//     {
//         std::cout << "[Init] Process 0: Generating Base FHE Keys..." << std::endl;
//         generator_key.keygen(0);
//         generator_key.save(
//             "client_secret.bin",
//             "server_eval.bin",
//             "server_annihilate.bin",
//             "server_pk.bin",
//             "server_pk_v2.bin");
//         std::cout << "[Init] Process 0: Keys saved successfully. Exiting." << std::endl;
//         return 0;
//     }

//     // ==========================================

//     // ==========================================
//     std::cout << "[Process " << process_id << "] Booting up. Loading keys..." << std::endl;
//     if (!is_exist(".", "client_secret.bin"))
//     {
//         std::cerr << "🚨 [Process " << process_id << "] Fatal: Keys not found! Run with ID 0 first." << std::endl;
//         return 1;
//     }
    
//     generator_key.load(
//         "server_eval.bin",
//         "server_annihilate.bin",
//         "server_pk.bin",
//         "server_pk_v2.bin");

//     map<string, trlwe_ptr> global_ct_map;
//     load_ciphertexts(global_ct_map, "data/ciphertexts_upload.txt", generator_key);

//     std::cout << "[Process " << process_id << "] Starting homomorphic computation..." << std::endl;
    

//     run(process_id, global_ct_map);

//     std::cout << "[Process " << process_id << "] Computation finished." << std::endl;
    

//     return 0;
// }
