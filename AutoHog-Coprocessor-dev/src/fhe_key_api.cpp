#include "fhe_key.h"
#include "fhe_key_api.h"

#include <cstdlib>
#include <cstring>
#include <limits>
#include <memory>
#include <stdexcept>

namespace {

constexpr int kSha256HexLength = 64;
constexpr int kMaxMetadataLength = 1 << 20;
constexpr uint64_t kMaxBulkElements = 4096;

void reset_out_len(int* out_len) noexcept {
    if (out_len != nullptr) *out_len = 0;
}

void require_buffer(const void* buffer, int len, const char* name,
                    bool allow_empty = false) {
    if (len < 0 || (!allow_empty && len == 0) || (len > 0 && buffer == nullptr)) {
        throw std::invalid_argument(std::string("invalid ") + name);
    }
}

Bytes copy_bytes(const uint8_t* data, int len, const char* name,
                 bool allow_empty = false) {
    require_buffer(data, len, name, allow_empty);
    return len == 0 ? Bytes{} : Bytes(data, data + len);
}

const uint8_t* copy_result(int* out_len, const Bytes& bytes) {
    if (out_len == nullptr) throw std::invalid_argument("out_len is null");
    if (bytes.empty() || bytes.size() > static_cast<size_t>(std::numeric_limits<int>::max())) {
        throw std::length_error("invalid serialized ciphertext length");
    }
    auto* result = static_cast<uint8_t*>(std::malloc(bytes.size()));
    if (result == nullptr) throw std::bad_alloc();
    std::memcpy(result, bytes.data(), bytes.size());
    *out_len = static_cast<int>(bytes.size());
    return result;
}

template<int VERSION>
std::unique_ptr<FheKey> create_key(const char* pk, int len) {
    require_buffer(pk, len, "public key");
    auto key = std::make_unique<FheKey>();
    if constexpr (VERSION == 1) key->load_pk(pk, len);
    else key->load_pk_v2(pk, len);
    return key;
}

std::unique_ptr<FheKey> create_key_with_zk(
    const char* pk, int pk_len, const char* crs, int crs_len,
    const char* expected_sha256, int expected_sha256_len) {
    require_buffer(crs, crs_len, "CRS");
    require_buffer(expected_sha256, expected_sha256_len, "expected CRS SHA-256");
    if (expected_sha256_len != kSha256HexLength) {
        throw std::invalid_argument("expected CRS SHA-256 must contain 64 hex characters");
    }
    auto key = create_key<1>(pk, pk_len);
    key->load_crs(crs, crs_len,
                  std::string(expected_sha256, expected_sha256 + expected_sha256_len));
    return key;
}

template<int VERSION>
const uint8_t* encrypt_integer_impl(int* out_len, FheKey* key,
                                    const uint8_t* in, int in_len) {
    if (key == nullptr) throw std::invalid_argument("key is null");
    constexpr size_t max_bytes =
        (VERSION == 1 ? lvl2param::n : lvl1param::n) / 8;
    if (in_len <= 0 || static_cast<size_t>(in_len) > max_bytes) {
        throw std::invalid_argument("plaintext length is outside the supported polynomial range");
    }
    Bytes plaintext = copy_bytes(in, in_len, "plaintext");
    Bytes serialized;
    if constexpr (VERSION == 1) {
        serialized = serialize_compressed_ciphertexts(key->encryptBatch(plaintext));
    } else {
        serialized = serialize_compressed_ciphertexts_v2(key->encryptBatchV2(plaintext));
    }
    return copy_result(out_len, serialized);
}

const uint8_t* encrypt_with_zk_impl(int* out_len, FheKey* key,
                                    const uint8_t* in, int in_len,
                                    const uint8_t* metadata, int metadata_len) {
    if (key == nullptr) throw std::invalid_argument("key is null");
    if (in_len <= 0 || static_cast<size_t>(in_len) > lvl2param::n / 8) {
        throw std::invalid_argument("plaintext length is outside the supported polynomial range");
    }
    if (metadata_len < 0 || metadata_len > kMaxMetadataLength) {
        throw std::invalid_argument("metadata length is invalid");
    }
    Bytes plaintext = copy_bytes(in, in_len, "plaintext");
    Bytes metadata_bytes = copy_bytes(metadata, metadata_len, "metadata", true);
    return copy_result(out_len, key->encryptWithZkBatch(plaintext, metadata_bytes));
}

uint64_t read_u64_be(const uint8_t*& cursor, size_t& remaining) {
    if (remaining < 8) throw std::invalid_argument("truncated bulk integer header");
    uint64_t value = 0;
    for (int i = 0; i < 8; ++i) value = (value << 8) | cursor[i];
    cursor += 8;
    remaining -= 8;
    return value;
}

void append_u64_be(Bytes& out, uint64_t value) {
    for (int shift = 56; shift >= 0; shift -= 8) {
        out.push_back(static_cast<uint8_t>((value >> shift) & 0xff));
    }
}

uint8_t value_type(uint64_t value_size) {
    switch (value_size) {
        case 1: return 1;
        case 2: return 2;
        case 4: return 3;
        case 8: return 4;
        case 16: return 5;
        case 32: return 6;
        default: throw std::invalid_argument("unsupported integer byte width");
    }
}

const uint8_t* encrypt_bulk_impl(int* out_len, FheKey* key,
                                 const uint8_t* in, int in_len,
                                 const uint8_t* metadata, int metadata_len) {
    if (key == nullptr) throw std::invalid_argument("key is null");
    require_buffer(in, in_len, "bulk plaintext");
    const uint8_t* cursor = in;
    size_t remaining = static_cast<size_t>(in_len);
    const uint64_t count = read_u64_be(cursor, remaining);
    if (count == 0 || count > kMaxBulkElements) {
        throw std::invalid_argument("bulk element count is invalid");
    }

    Bytes output;
    append_u64_be(output, count);
    for (uint64_t i = 0; i < count; ++i) {
        const uint64_t value_size = read_u64_be(cursor, remaining);
        const uint8_t type = value_type(value_size);
        if (value_size > remaining || value_size > static_cast<uint64_t>(std::numeric_limits<int>::max())) {
            throw std::invalid_argument("truncated bulk integer payload");
        }

        int ct_len = 0;
        std::unique_ptr<uint8_t, decltype(&std::free)> ct(
            const_cast<uint8_t*>(encrypt_with_zk_impl(
                &ct_len, key, cursor, static_cast<int>(value_size),
                metadata, metadata_len)), &std::free);
        cursor += value_size;
        remaining -= static_cast<size_t>(value_size);

        const size_t extra = 1 + 8 + static_cast<size_t>(ct_len);
        if (extra > static_cast<size_t>(std::numeric_limits<int>::max()) - output.size()) {
            throw std::length_error("bulk ciphertext is too large");
        }
        output.push_back(type);
        append_u64_be(output, static_cast<uint64_t>(ct_len));
        output.insert(output.end(), ct.get(), ct.get() + ct_len);
    }
    if (remaining != 0) throw std::invalid_argument("trailing bytes in bulk plaintext");
    return copy_result(out_len, output);
}

} // namespace

extern "C" {

void* create_fhe_key(const char* pk, int len) {
    try { return create_key<1>(pk, len).release(); }
    catch (...) { return nullptr; }
}

void* create_fhe_key_v2(const char* pk, int len) {
    try { return create_key<2>(pk, len).release(); }
    catch (...) { return nullptr; }
}

void* create_fhe_key_with_zk(const char*, int, const char*, int) {
    return nullptr;
}

void* create_fhe_key_with_zk_sha256(const char* pk, int pk_len,
                                     const char* crs, int crs_len,
                                     const char* expected_sha256,
                                     int expected_sha256_len) {
    try {
        return create_key_with_zk(pk, pk_len, crs, crs_len,
                                  expected_sha256, expected_sha256_len).release();
    } catch (...) { return nullptr; }
}

void free_fhe_key(const void* key) {
    try { delete static_cast<const FheKey*>(key); }
    catch (...) {}
}

const uint8_t* encrypt_integer(int* out_len, void* key,
                               const uint8_t* in, int in_len) {
    reset_out_len(out_len);
    try { return encrypt_integer_impl<1>(out_len, static_cast<FheKey*>(key), in, in_len); }
    catch (...) { return nullptr; }
}

const uint8_t* encrypt_integer_v2(int* out_len, void* key,
                                  const uint8_t* in, int in_len) {
    reset_out_len(out_len);
    try { return encrypt_integer_impl<2>(out_len, static_cast<FheKey*>(key), in, in_len); }
    catch (...) { return nullptr; }
}

const uint8_t* encrypt_integer_ex(int* out_len, const char* pk, int len,
                                  const uint8_t* in, int in_len) {
    reset_out_len(out_len);
    try {
        auto key = create_key<1>(pk, len);
        return encrypt_integer_impl<1>(out_len, key.get(), in, in_len);
    } catch (...) { return nullptr; }
}

const uint8_t* encrypt_integer_v2_ex(int* out_len, const char* pk, int len,
                                     const uint8_t* in, int in_len) {
    reset_out_len(out_len);
    try {
        auto key = create_key<2>(pk, len);
        return encrypt_integer_impl<2>(out_len, key.get(), in, in_len);
    } catch (...) { return nullptr; }
}

const uint8_t* encrypt_bool(int* out_len, void* key, bool value) {
    const uint8_t byte = value ? 1 : 0;
    return encrypt_integer(out_len, key, &byte, 1);
}

const uint8_t* encrypt_bool_v2(int* out_len, void* key, bool value) {
    const uint8_t byte = value ? 1 : 0;
    return encrypt_integer_v2(out_len, key, &byte, 1);
}

const uint8_t* encrypt_bool_ex(int* out_len, const char* pk, int len, bool value) {
    const uint8_t byte = value ? 1 : 0;
    return encrypt_integer_ex(out_len, pk, len, &byte, 1);
}

const uint8_t* encrypt_bool_v2_ex(int* out_len, const char* pk, int len, bool value) {
    const uint8_t byte = value ? 1 : 0;
    return encrypt_integer_v2_ex(out_len, pk, len, &byte, 1);
}

const uint8_t* encrypt_integer_with_zk(int* out_len, void* key,
                                       const uint8_t* in, int in_len,
                                       const uint8_t* metadata, int metadata_len) {
    reset_out_len(out_len);
    try {
        return encrypt_with_zk_impl(out_len, static_cast<FheKey*>(key), in, in_len,
                                    metadata, metadata_len);
    } catch (...) { return nullptr; }
}

const uint8_t* encrypt_bool_with_zk(int* out_len, void* key, bool value,
                                    const uint8_t* metadata, int metadata_len) {
    const uint8_t byte = value ? 1 : 0;
    return encrypt_integer_with_zk(out_len, key, &byte, 1, metadata, metadata_len);
}

const uint8_t* encrypt_multi_integers_with_zk(int* out_len, void* key,
                                              const uint8_t* in, int in_len,
                                              const uint8_t* metadata, int metadata_len) {
    reset_out_len(out_len);
    try {
        return encrypt_bulk_impl(out_len, static_cast<FheKey*>(key), in, in_len,
                                 metadata, metadata_len);
    } catch (...) { return nullptr; }
}

const uint8_t* encrypt_multi_integers_with_zk_ex(
    int* out_len, const char*, int, const char*, int,
    const uint8_t*, int, const uint8_t*, int) {
    reset_out_len(out_len);
    return nullptr;
}

const uint8_t* encrypt_multi_integers_with_zk_ex_sha256(
    int* out_len, const char* pk, int pk_len,
    const char* crs, int crs_len,
    const char* expected_sha256, int expected_sha256_len,
    const uint8_t* in, int in_len,
    const uint8_t* metadata, int metadata_len) {
    reset_out_len(out_len);
    try {
        auto key = create_key_with_zk(pk, pk_len, crs, crs_len,
                                      expected_sha256, expected_sha256_len);
        return encrypt_bulk_impl(out_len, key.get(), in, in_len, metadata, metadata_len);
    } catch (...) { return nullptr; }
}

void free_data(const uint8_t* data) {
    std::free(const_cast<uint8_t*>(data));
}

} // extern "C"
