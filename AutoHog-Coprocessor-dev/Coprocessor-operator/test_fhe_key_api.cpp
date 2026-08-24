#include "fhe_key.h"
#include "fhe_key_api.h"

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

    FheKey* api_key = nullptr;
    {
        FILE *fp = fopen("server_pk.bin", "rb");
        if (!fp) {
            cerr<< "open file error" << endl;
            exit(1);
        }
        fseek(fp, 0, SEEK_END);
        int len = ftell(fp);
        fseek(fp, 0, SEEK_SET);
        char* buf = new char[len];
        fread(buf, 1, len, fp);
        api_key = (FheKey*)create_fhe_key(buf, len);
        fclose(fp);
    }

    vector<uint8_t> lhs(32, 0);
    lhs[31] = 1;
    vector<uint8_t> rhs(32, 0xff);
    // rhs[0] = 0x7f;
    auto encrypt_fn = [&api_key](const vector<uint8_t>& data) {
        int out_len;
        const uint8_t* buf = encrypt_integer(&out_len, api_key, data.data(), data.size());
        vector<uint8_t> vec(buf, buf + out_len);
        auto ct = deserialize_compressed_ciphertexts(vec);
        free_data(buf);
        return ct;
    };

    auto lhs_ct = decompress_ciphertexts(encrypt_fn(lhs));
    auto rhs_ct = decompress_ciphertexts(encrypt_fn(rhs));

    free_fhe_key(api_key);

    auto res_ct = key.lt(lhs_ct, rhs_ct);
    auto bit_ct = to_tlwe_ptr(res_ct);
    auto bit = key.decrypt(bit_ct);
    cout << "lt result: " << (bit? "1": "0") << endl;

    return 0;
}
