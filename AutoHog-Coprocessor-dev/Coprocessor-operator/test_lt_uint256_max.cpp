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
        key.keygen(0);
        key.save("client_secret.bin", "server_eval.bin", "server_ahk.bin", "server_pk.bin", "server_pk_v2.bin");
    }
    {
        key.load("server_eval.bin", "server_ahk.bin", "server_pk.bin", "server_pk_v2.bin");
	key.load_sk("client_secret.bin");
    }

    vector<uint8_t> lhs(32, 0);
    lhs[31] = 1;
    vector<uint8_t> rhs(32, 0xff);
    // rhs[0] = 0x7f;

    auto lhs_ct = decompress_ciphertexts(key.trivialEncryptBatch(lhs));
    auto rhs_ct = decompress_ciphertexts(key.trivialEncryptBatch(rhs));

    auto res_ct = key.lt(lhs_ct, rhs_ct);
    auto bit_ct = to_tlwe_ptr(res_ct);
    auto bit = key.decrypt(bit_ct);
    cout << "lt result: " << (bit? "1": "0") << endl;

    return 0;
}
