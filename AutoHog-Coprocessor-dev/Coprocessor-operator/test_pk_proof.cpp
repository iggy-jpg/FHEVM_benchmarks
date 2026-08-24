#include "fhe_key.h"
#include <cstdlib>

bool is_exist(string dirname, string filename) {
    std::ifstream file(dirname + "/" + filename, std::ios::binary);
    if (file.is_open()) {
        file.close();
        return true;
    }
    return false;
}

std::string trusted_crs_sha256() {
    const char* value = std::getenv("AUTOHOG_CRS_SHA256");
    if (value == nullptr || std::strlen(value) != 64) {
        throw std::runtime_error("AUTOHOG_CRS_SHA256 must contain the pinned 64-hex CRS digest");
    }
    return value;
}

void load_crs_file(FheKey& key, const char* filename, const std::string& digest) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        throw std::runtime_error(std::string("missing CRS file: ") + filename);
    }
    std::vector<char> crs(
        (std::istreambuf_iterator<char>(file)),
        std::istreambuf_iterator<char>()
    );
    if (crs.empty()) {
        throw std::runtime_error(std::string("empty CRS file: ") + filename);
    }
    key.load_crs(crs.data(), static_cast<int>(crs.size()), digest);
}

int main(int argc, const char* argv[]) {
    const std::string crs_digest = trusted_crs_sha256();
    FheKey key;
    if (!is_exist(".", "client_secret.bin")) {
        key.keygen(4);
        load_crs_file(key, "crs.bin", crs_digest);
        key.save("client_secret.bin", "server_eval.bin", "server_ahk.bin", "server_pk.bin", "server_pk_v2.bin", "crs.bin");
    }
    {
        key.load("server_eval.bin", "server_ahk.bin", "server_pk.bin", "server_pk_v2.bin", "crs.bin", crs_digest);
	    key.load_sk("client_secret.bin");
    }
    const std::vector<uint8_t> metadata = {'M', 'E', 'T', 'A', 'D', 'A', 'T', 'A'};

    vector<uint8_t> lhs(32, 0);
    lhs[31] = 1;
    vector<uint8_t> rhs(32, 0);
    rhs[31] = 2;
    uint8_t ctl = 0;

    auto test_add = [&](const trlwe_ptr& lhs0, const trlwe_ptr& rhs0) {
        auto lhs_ct = decompress_ciphertexts(lhs0);
        auto rhs_ct = decompress_ciphertexts(rhs0);
    
        auto res_ct = key.add(lhs_ct, rhs_ct);

        auto print = [&](tlwe_vec_ptr& ct) {
            auto bytes = key.decryptBatch(ct);
            std::for_each(bytes.begin(), bytes.end(), [](uint8_t byte) {
                printf("%02x ", byte);
            });
            printf("\n");
        };

        print(lhs_ct);
        print(rhs_ct);
        print(res_ct);
    };
    auto test_select = [&](const trlwe_ptr& lhs0, const trlwe_ptr& rhs0, const trlwe_ptr& ctl) {
        auto lhs_ct = decompress_ciphertexts(lhs0);
        auto rhs_ct = decompress_ciphertexts(rhs0);
        auto ctl_ct = decompress_ciphertexts(ctl);
    
        auto res_ct = key.select(lhs_ct, rhs_ct, ctl_ct);

        auto print = [&](tlwe_vec_ptr& ct) {
            auto bytes = key.decryptBatch(ct);
            std::for_each(bytes.begin(), bytes.end(), [](uint8_t byte) {
                printf("%02x ", byte);
            });
            printf("\n");
        };

        print(lhs_ct);
        print(rhs_ct);
        print(ctl_ct);
        print(res_ct);
    };
        

    auto test_encrypt_batch = [&](const vector<uint8_t>& data) {
        return key.encryptWithZkBatch(data, metadata);
    };

    auto test_encrypt = [&](uint8_t data) {
        return key.encryptWithZk(data, metadata);
    };

    auto test_verify = [&](const vector<uint8_t>& data) {
        trlwe_ptr ct;
        bool ok = key.verifyZkProof(ct, data, metadata);
        if (!ok) {
            throw std::runtime_error("verify failed");
        }
        return std::move(ct);
    };

    auto lhs0 = test_encrypt_batch(lhs);
    auto rhs0 = test_encrypt_batch(rhs);
    auto ctl0 = test_encrypt(ctl);
    printf("%d\n", lhs0.size());

    auto lhs1 = test_verify(lhs0);
    auto rhs1 = test_verify(rhs0);
    auto ctl1 = test_verify(ctl0);

    test_add(lhs1, rhs1);
    printf("\n");
    test_select(lhs1, rhs1, ctl1);

    return 0;
}
