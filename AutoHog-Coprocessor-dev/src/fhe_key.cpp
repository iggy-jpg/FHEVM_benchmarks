#include "fhe_key.h"
#include "zk_backend.h"

void load_pk(PKE::SeededRlwePkPair<lvl2param> *spk, const char* pk, int len){
    if (spk == nullptr || pk == nullptr || len < 0) {
        throw std::invalid_argument("invalid public key buffer");
    }
    size_t seed_len = sizeof(spk->seed);
    size_t b_len = sizeof(typename lvl2param::T) * lvl2param::n;
    size_t ct_zero_len = sizeof(CT_ZERO);
    size_t ct_one_len = sizeof(CT_ONE);
    if (len != seed_len + b_len + ct_zero_len + ct_one_len) {
        throw std::runtime_error("error pk len");
    }
    memcpy(&spk->seed, pk, seed_len);
    pk += seed_len;
    memcpy((char*)spk->b.data(), pk, b_len);
    pk += b_len;
    // Serialized constants are retained only for wire-format compatibility.
    // They are caller-controlled and must never overwrite process-global state.
    pk += ct_zero_len;
    (void)pk;
}
void load_pk_v2(PKE::SeededRlwePkPair<lvl1param> *spk, const char* pk, int len){
    if (spk == nullptr || pk == nullptr || len < 0) {
        throw std::invalid_argument("invalid public key buffer");
    }
    size_t seed_len = sizeof(spk->seed);
    size_t b_len = sizeof(typename lvl1param::T) * lvl1param::n;
    if (len != seed_len + b_len) {
        throw std::runtime_error("error pk len");
    }
    memcpy(&spk->seed, pk, seed_len);
    pk += seed_len;
    memcpy((char*)spk->b.data(), pk, b_len);
    pk += b_len;
}

void load_crs(::ZKProof::PublicParams<crs_lvl>* params, const char* crs, int len) {
    (void)params;
    (void)crs;
    (void)len;
    throw std::invalid_argument(
        "Loading CRS bytes requires an independently trusted SHA-256 digest");
}

void load_crs(::ZKProof::PublicParams<crs_lvl>* params, const char* crs, int len,
              const std::string& expected_sha256) {
    if (params == nullptr) throw std::invalid_argument("CRS destination is null");
    *params = ZKProof::load_public_params_from_memory_by_expected_hash<crs_lvl>(
        crs, len, expected_sha256);
}

Bytes serialize_ciphertexts(const tlwe_vec_ptr& d) {
    require_valid_tlwe_vector(d, "TLWE vector to serialize");
    Bytes tmp(d->data.size() * sizeof(tlwe_t));
    const uint8_t* p = reinterpret_cast<const uint8_t*>(d->data.data());
    memcpy(tmp.data(), p, tmp.size());
    return tmp;
}

template<typename param>
Bytes serialize_compressed_ciphertexts_generic(const std::unique_ptr<TRLWEWrapperGeneric<param>>& d) {
    require_valid_trlwe_wrapper<param>(d, "compressed ciphertext to serialize");
    size_t element_size = sizeof(typename param::T) * param::n;
    size_t num_elements = param::k + 1;
    vector<uint8_t> vec(4 + num_elements * element_size); // the first 4 bytes is W
    uint32_t W = d->W;
    memcpy(&vec[0], &W, sizeof(W));

    for (size_t i = 0; i < num_elements; i++) {
        memcpy((&vec[4]) + i * element_size, d->data[i].data(), element_size);
    }
    return vec;
}

Bytes serialize_compressed_ciphertexts(const trlwe_ptr& d) {
    return serialize_compressed_ciphertexts_generic<lvl2param>(d);
}
Bytes serialize_compressed_ciphertexts_v2(const trlwe_v2_ptr& d) {
    return serialize_compressed_ciphertexts_generic<lvl1param>(d);
}

static void save_into_file(const string& file, const Bytes& bytes) {
    if (bytes.empty()) throw std::invalid_argument("cannot save an empty ciphertext");
    std::ofstream ofs(file, std::ios::binary);
    if (!ofs) throw std::runtime_error("cannot open ciphertext file for writing");
    ofs.write(reinterpret_cast<const char*>(bytes.data()),
              static_cast<std::streamsize>(bytes.size()));
    ofs.flush();
    if (!ofs) throw std::runtime_error("failed writing ciphertext file");
    ofs.close();
}
void save_ciphertexts(const string& file, const tlwe_vec_ptr& d) {
    auto bytes = serialize_ciphertexts(d);
    bytes.insert(bytes.begin(), 1);
    save_into_file(file, bytes);
}
void save_compressed_ciphertexts(const string& file, const trlwe_ptr& d) {
    auto bytes = serialize_compressed_ciphertexts(d);
    bytes.insert(bytes.begin(), 2);
    save_into_file(file, bytes);
}
void save_compressed_ciphertexts_v2(const string& file, const trlwe_v2_ptr& d) {
    auto bytes = serialize_compressed_ciphertexts_v2(d);
    bytes.insert(bytes.begin(), 3);
    save_into_file(file, bytes);
}
tlwe_vec_ptr deserialize_ciphertexts(const Bytes& bytes) {
    if (bytes.empty() || bytes.size() % sizeof(tlwe_t) != 0) {
        throw std::invalid_argument("invalid TLWE ciphertext byte length");
    }
    size_t num = bytes.size() / sizeof(tlwe_t);
    if (num > lvl2param::n) {
        throw std::length_error("too many TLWE ciphertexts");
    }
    tlwe_vec_t ciphertexts(num);
    uint8_t* p = reinterpret_cast<uint8_t*>(ciphertexts.data());
    memcpy(p, bytes.data(), bytes.size());
    return tlwe_vec_ptr(new tlwe_vec(ciphertexts));
}

template<typename param>
std::unique_ptr<TRLWEWrapperGeneric<param>> deserialize_compressed_ciphertexts_generic(const Bytes& bytes) {
    size_t element_size = sizeof(typename param::T) * param::n;
    size_t num_elements = param::k + 1;
    const size_t expected_size = sizeof(uint32_t) + num_elements * element_size;
    if (bytes.size() != expected_size) {
        throw std::invalid_argument("invalid compressed ciphertext byte length");
    }
    TRLWE<param> result;
    uint32_t serialized_w = 0;
    std::memcpy(&serialized_w, bytes.data(), sizeof(serialized_w));
    if (serialized_w == 0 || serialized_w > param::n) {
        throw std::invalid_argument("invalid compressed ciphertext bit width");
    }
    size_t W = serialized_w;

    for (size_t i = 0; i < num_elements; i++) {
        memcpy((uint8_t*)result[i].data(), (&bytes[4]) + i * element_size, element_size);
    }
    return std::unique_ptr<TRLWEWrapperGeneric<param>>(new TRLWEWrapperGeneric<param>(result, W));
}

trlwe_ptr deserialize_compressed_ciphertexts(const Bytes& bytes) {
    return deserialize_compressed_ciphertexts_generic<lvl2param>(bytes);
}
trlwe_v2_ptr deserialize_compressed_ciphertexts_v2(const Bytes& bytes) {
    return deserialize_compressed_ciphertexts_generic<lvl1param>(bytes);
}
Bytes load_from_file(const string& file, size_t max_size) {
    std::ifstream ifs(file, std::ios::binary);
    if (!ifs) throw std::runtime_error("cannot open ciphertext file");
    ifs.seekg(0, std::ios::end);
    const std::streamoff end = ifs.tellg();
    if (end <= 0 || static_cast<uint64_t>(end) > max_size) {
        throw std::length_error("ciphertext file size is invalid");
    }
    size_t size = static_cast<size_t>(end);
    ifs.seekg(0, std::ios::beg);

    Bytes bytes(size);
    ifs.read(reinterpret_cast<char*>(bytes.data()), static_cast<std::streamsize>(bytes.size()));
    if (!ifs || ifs.gcount() != static_cast<std::streamsize>(bytes.size())) {
        throw std::runtime_error("truncated ciphertext file");
    }
    ifs.close();
    return bytes;
}

tlwe_vec_ptr load_ciphertexts(const string& file) {
    const size_t max_size = 1 + lvl2param::n * sizeof(tlwe_t);
    auto bytes = load_from_file(file, max_size);
    if (bytes[0] != 1) {
        throw std::runtime_error("ciphertext format error");
    }
    bytes.erase(bytes.begin());
    return deserialize_ciphertexts(bytes);
}
trlwe_ptr load_compressed_ciphertexts(const string& file) {
    const size_t exact_size = 1 + sizeof(uint32_t) +
        (lvl2param::k + 1) * sizeof(lvl2param::T) * lvl2param::n;
    auto bytes = load_from_file(file, exact_size);
    if (bytes.size() != exact_size) throw std::invalid_argument("ciphertext file size mismatch");
    if (bytes[0] != 2) {
        throw std::runtime_error("ciphertext format error");
    }
    bytes.erase(bytes.begin());
    return deserialize_compressed_ciphertexts(bytes);
}
trlwe_v2_ptr load_compressed_ciphertexts_v2(const string& file) {
    const size_t exact_size = 1 + sizeof(uint32_t) +
        (lvl1param::k + 1) * sizeof(lvl1param::T) * lvl1param::n;
    auto bytes = load_from_file(file, exact_size);
    if (bytes.size() != exact_size) throw std::invalid_argument("ciphertext file size mismatch");
    if (bytes[0] != 3) {
        throw std::runtime_error("ciphertext format error");
    }
    bytes.erase(bytes.begin());
    return deserialize_compressed_ciphertexts_v2(bytes);
}

template<typename P>
Bytes serialize_ciphertext_with_proof_generic(const trlwe_ptr& ct, const ZKProof::Proof& proof) {
    require_valid_trlwe_wrapper<lvl2param>(ct, "ciphertext with proof");
    stringstream ss;
    // Magic
    const char magic[4] = {'Z', 'K', 'S', '1'};
    ss.write(magic, 4);

    // Ciphertext
    for (size_t i = 0; i < P::k + 1; ++i) {
        ss.write(reinterpret_cast<const char*>(ct->data[i].data()),
                  static_cast<std::streamsize>(sizeof(typename P::T) * P::n));
    }

    uint64_t w = ct->W;
    ss.write(reinterpret_cast<const char*>(&w), sizeof(w));

    // Proof
    auto proof_bytes = proof.serialize();
    uint64_t proof_size = proof_bytes.size();
    ss.write(reinterpret_cast<const char*>(&proof_size), sizeof(proof_size));
    ss.write(reinterpret_cast<const char*>(proof_bytes.data()),
              static_cast<std::streamsize>(proof_size));

    string buffer = ss.str();
    std::vector<uint8_t> vec(buffer.begin(), buffer.end());
    return vec;
}


Bytes serialize_ciphertext_with_proof(const trlwe_ptr& ciphertext, const ZKProof::Proof& proof) {
    return serialize_ciphertext_with_proof_generic<crs_lvl>(ciphertext, proof);
}
template<typename P>
bool deserialize_ciphertext_with_proof_generic(trlwe_ptr& ciphertext, ZKProof::Proof& proof, const Bytes& data) {
    try {
        constexpr uint64_t kMaxProofBytes = 1ULL << 20;
        const size_t ciphertext_bytes = (P::k + 1) * sizeof(typename P::T) * P::n;
        const size_t fixed_bytes = 4 + ciphertext_bytes + sizeof(uint64_t) + sizeof(uint64_t);
        if (data.size() < fixed_bytes) return false;
        size_t offset = 0;
        trlwe_t ct;

        if (std::memcmp(data.data(), "ZKS1", 4) != 0) {
            std::cerr << "[loadSealedCiphertext] Invalid file format (bad magic)\n";
            return false;
        }
        offset += 4;

        for (size_t i = 0; i < P::k + 1; ++i) {
            const size_t component_bytes = sizeof(typename P::T) * P::n;
            std::memcpy(ct[i].data(), data.data() + offset, component_bytes);
            offset += component_bytes;
        }
        uint64_t w = 0;
        std::memcpy(&w, data.data() + offset, sizeof(w));
        offset += sizeof(w);
        if (w == 0 || w > P::n) return false;

        uint64_t proof_size = 0;
        std::memcpy(&proof_size, data.data() + offset, sizeof(proof_size));
        offset += sizeof(proof_size);
        if (proof_size == 0 || proof_size > kMaxProofBytes ||
            proof_size != data.size() - offset) return false;

        std::vector<uint8_t> proof_bytes(data.begin() + static_cast<std::ptrdiff_t>(offset),
                                         data.end());
        ZKProof::Proof parsed_proof = ZKProof::Proof::deserialize(proof_bytes);
        trlwe_ptr parsed_ciphertext(new trlwe(ct, w));
        ciphertext = std::move(parsed_ciphertext);
        proof = std::move(parsed_proof);

        return true;
    } catch (const std::exception& e) {
        std::cerr << "[loadSealedCiphertext] Exception: " << e.what() << "\n";
        return false;
    }
}

bool deserialize_ciphertext_with_proof(trlwe_ptr& ct, ZKProof::Proof& proof, const Bytes& data) {
    return deserialize_ciphertext_with_proof_generic<crs_lvl>(ct, proof, data);
}
