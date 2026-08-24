#include "fhe_key.h"

bool is_exist(string dirname, string filename) {
    std::ifstream file(dirname + "/" + filename, std::ios::binary);
    if (file.is_open()) {
        file.close();
        return true;
    }
    return false;
}

int main(int argc, const char* argv[]) {
    FheKey key;
    if (!is_exist(".", "client_secret.bin")) {
        key.keygen(4);
        key.save("client_secret.bin", "server_eval.bin", "server_ahk.bin", "server_pk.bin", "server_pk_v2.bin");
    }
    {
        key.load("server_eval.bin", "server_ahk.bin", "server_pk.bin", "server_pk_v2.bin");
	    key.load_sk("client_secret.bin");
    }

    vector<uint8_t> lhs(32, 0);
    lhs[31] = 1;
    vector<uint8_t> rhs(32, 0);
    rhs[31] = 2;

    auto test_pk = [&](bool use_pk_v2) {
        std::function<trlwe_ptr (vector<uint8_t>)> encrypt_fn;
        if (use_pk_v2) {
            encrypt_fn = [&](vector<uint8_t> data) { return key.decompressV3(key.encryptBatchV2(data)); };
        }
        else {
            encrypt_fn = [&](vector<uint8_t> data) {return key.encryptBatch(data);};
        }
        auto lhs_ct = decompress_ciphertexts(encrypt_fn(lhs));
        auto rhs_ct = decompress_ciphertexts(encrypt_fn(rhs));
    
        auto res_ct = key.add(lhs_ct, rhs_ct);
        auto bytes = key.decryptBatch(res_ct);
        printf("use pk v2: %d ", use_pk_v2);
        std::for_each(bytes.begin(), bytes.end(), [](uint8_t byte) {
            printf("%02x ", byte);
        });
        printf("\n");
    };

    test_pk(false);
    test_pk(true);

    return 0;
}
