#pragma once
#include <cstdint>
#include <string>
#include <cstring>

extern "C" {
void* create_fhe_key(const char* pk, int len);
void* create_fhe_key_v2(const char* pk, int len);
// Legacy unauthenticated entry point: retained for ABI compatibility and fails closed.
void* create_fhe_key_with_zk(const char* pk, int pk_len, const char* crs, int crs_len);
void* create_fhe_key_with_zk_sha256(const char* pk, int pk_len,
                                     const char* crs, int crs_len,
                                     const char* expected_sha256,
                                     int expected_sha256_len);

void free_fhe_key(const void* key);

const uint8_t* encrypt_integer(int* out_len, void* key, const uint8_t* in, int in_len);
const uint8_t* encrypt_integer_ex(int* out_len, const char* pk, int len, const uint8_t* in, int in_len);
const uint8_t* encrypt_integer_v2(int* out_len, void* key, const uint8_t* in, int in_len);
const uint8_t* encrypt_integer_v2_ex(int* out_len, const char* pk, int len, const uint8_t* in, int in_len);
const uint8_t* encrypt_integer_with_zk(int* out_len, void* key, const uint8_t* in, int in_len, const uint8_t* metadata, int metadata_len);
const uint8_t* encrypt_multi_integers_with_zk(int* out_len, void* key, const uint8_t* in, int in_len, const uint8_t* metadata, int metadata_len);
// Legacy unauthenticated entry point: retained for ABI compatibility and fails closed.
const uint8_t* encrypt_multi_integers_with_zk_ex(int* out_len, const char* pk, int pk_len, const char* crs, int crs_len, const uint8_t* in, int in_len, const uint8_t* metadata, int metadata_len);
const uint8_t* encrypt_multi_integers_with_zk_ex_sha256(
    int* out_len, const char* pk, int pk_len,
    const char* crs, int crs_len,
    const char* expected_sha256, int expected_sha256_len,
    const uint8_t* in, int in_len,
    const uint8_t* metadata, int metadata_len);

const uint8_t* encrypt_bool(int* out_len, void* key, bool b); 

const uint8_t* encrypt_bool_ex(int* out_len, const char* pk, int len, bool b); 
const uint8_t* encrypt_bool_v2(int* out_len, void* key, bool b); 

const uint8_t* encrypt_bool_v2_ex(int* out_len, const char* pk, int len, bool b); 
const uint8_t* encrypt_bool_with_zk(int* out_len, void* key, bool b, const uint8_t* metadata, int metadata_len); 

void free_data(const uint8_t* data);
}
