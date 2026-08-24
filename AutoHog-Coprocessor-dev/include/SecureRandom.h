#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>
#include <type_traits>

#if !defined(__EMSCRIPTEN__)
#include <openssl/rand.h>
#endif
#include <openssl/evp.h>

namespace autohog {

inline void secure_erase(void* ptr, std::size_t len)
{
    volatile std::uint8_t* bytes = static_cast<volatile std::uint8_t*>(ptr);
    while (len-- != 0) {
        *bytes++ = 0;
    }
}

// Fill the requested buffer from a cryptographic entropy source. Failure is
// fatal: cryptographic callers must never silently fall back to a weak PRNG.
inline void secure_random_bytes(std::uint8_t* out, std::size_t len)
{
    if (out == nullptr && len != 0) {
        throw std::invalid_argument("secure_random_bytes received a null buffer");
    }

#if defined(__EMSCRIPTEN__)
    // libc++'s random_device is backed by crypto.getRandomValues in Emscripten.
    std::random_device rd;
    while (len != 0) {
        const std::uint32_t word = static_cast<std::uint32_t>(rd());
        const std::size_t take = std::min(len, sizeof(word));
        std::memcpy(out, &word, take);
        out += take;
        len -= take;
    }
#else
    while (len != 0) {
        const std::size_t max_chunk =
            static_cast<std::size_t>((std::numeric_limits<int>::max)());
        const int chunk = static_cast<int>(std::min(len, max_chunk));
        if (RAND_bytes(out, chunk) != 1) {
            throw std::runtime_error("RAND_bytes failed while generating secure randomness");
        }
        out += chunk;
        len -= static_cast<std::size_t>(chunk);
    }
#endif
}

// Public seeds are serialized with public keys. 16 bytes provides the
// required 128 bits of entropy without conflating the seed with secret noise.
using PublicSeed = std::array<std::uint8_t, 16>;
static_assert(sizeof(PublicSeed) == 16, "PublicSeed must contain 128 bits");

inline PublicSeed secure_public_seed()
{
    PublicSeed seed{};
    secure_random_bytes(seed.data(), seed.size());
    return seed;
}

inline void store_u32_be(std::uint8_t* out, std::uint32_t value)
{
    out[0] = static_cast<std::uint8_t>(value >> 24);
    out[1] = static_cast<std::uint8_t>(value >> 16);
    out[2] = static_cast<std::uint8_t>(value >> 8);
    out[3] = static_cast<std::uint8_t>(value);
}

inline void store_u64_be(std::uint8_t* out, std::uint64_t value)
{
    for (int i = 7; i >= 0; --i) {
        out[i] = static_cast<std::uint8_t>(value);
        value >>= 8;
    }
}

// Deterministic public-mask expansion following the same construction pattern
// as TFHE-rs: a 128-bit compression seed is used as an AES-128 key and every
// object consumes a disjoint region of one counter-mode byte stream.
//
// Counter block layout (fixed protocol encoding):
//   bytes 0..3   = "AHG1" (format/version marker)
//   bytes 4..7   = big-endian domain identifier
//   bytes 8..15  = big-endian AES block index
//
// AES-ECB(counter_block) is the standard implementation of the CTR keystream;
// no plaintext is encrypted here. Public values are decoded little-endian so
// reconstruction is independent of host endianness and C++ standard library.
class PublicSeedStream {
public:
    PublicSeedStream(const PublicSeed& seed,
                     std::uint32_t domain,
                     std::uint64_t byte_offset = 0)
        : context_(EVP_CIPHER_CTX_new(), &EVP_CIPHER_CTX_free),
          domain_(domain),
          next_block_index_(byte_offset / 16u)
    {
        if (!context_) {
            throw std::runtime_error("EVP_CIPHER_CTX_new failed for public seed stream");
        }
        if (EVP_EncryptInit_ex(context_.get(), EVP_aes_128_ecb(), nullptr,
                               seed.data(), nullptr) != 1 ||
            EVP_CIPHER_CTX_set_padding(context_.get(), 0) != 1) {
            throw std::runtime_error("AES-128 initialization failed for public seed stream");
        }

        refill();
        offset_ = static_cast<std::size_t>(byte_offset % block_.size());
    }

    PublicSeedStream(const PublicSeedStream&) = delete;
    PublicSeedStream& operator=(const PublicSeedStream&) = delete;
    PublicSeedStream(PublicSeedStream&&) noexcept = default;
    PublicSeedStream& operator=(PublicSeedStream&&) noexcept = default;

    void fill_bytes(std::uint8_t* out, std::size_t len)
    {
        if (out == nullptr && len != 0) {
            throw std::invalid_argument("PublicSeedStream::fill_bytes received a null buffer");
        }

        while (len != 0) {
            if (offset_ == block_.size()) {
                refill();
            }
            const std::size_t take = std::min(len, block_.size() - offset_);
            std::memcpy(out, block_.data() + offset_, take);
            out += take;
            len -= take;
            offset_ += take;
        }
    }

    template<class UInt>
    UInt next_uniform()
    {
        static_assert(std::is_integral<UInt>::value && std::is_unsigned<UInt>::value,
                      "public coefficients must be unsigned integral values");
        static_assert(sizeof(UInt) <= sizeof(std::uint64_t),
                      "public coefficient wider than 64 bits is unsupported");

        std::array<std::uint8_t, sizeof(UInt)> bytes{};
        fill_bytes(bytes.data(), bytes.size());
        UInt value = 0;
        for (std::size_t i = 0; i < bytes.size(); ++i) {
            value |= static_cast<UInt>(bytes[i]) << (8 * i);
        }
        return value;
    }

private:
    using EvpContextPtr = std::unique_ptr<EVP_CIPHER_CTX, decltype(&EVP_CIPHER_CTX_free)>;

    void refill()
    {
        std::array<std::uint8_t, 16> counter{
            static_cast<std::uint8_t>('A'), static_cast<std::uint8_t>('H'),
            static_cast<std::uint8_t>('G'), static_cast<std::uint8_t>('1')
        };
        store_u32_be(counter.data() + 4, domain_);
        store_u64_be(counter.data() + 8, next_block_index_);

        int produced = 0;
        if (EVP_EncryptUpdate(context_.get(), block_.data(), &produced,
                              counter.data(), static_cast<int>(counter.size())) != 1 ||
            produced != static_cast<int>(block_.size())) {
            throw std::runtime_error("AES-CTR block expansion failed for public seed stream");
        }
        if (next_block_index_ == (std::numeric_limits<std::uint64_t>::max)()) {
            throw std::overflow_error("public seed stream exhausted");
        }
        ++next_block_index_;
        offset_ = 0;
    }

    EvpContextPtr context_;
    std::uint32_t domain_;
    std::uint64_t next_block_index_;
    std::array<std::uint8_t, 16> block_{};
    std::size_t offset_ = block_.size();
};

template<class UInt>
inline std::uint64_t public_object_byte_offset(std::uint64_t object_index,
                                               std::size_t elements_per_object)
{
    static_assert(std::is_integral<UInt>::value && std::is_unsigned<UInt>::value,
                  "public coefficients must be unsigned integral values");
    const std::uint64_t bytes_per_element = sizeof(UInt);
    if (elements_per_object >
        (std::numeric_limits<std::uint64_t>::max)() / bytes_per_element) {
        throw std::overflow_error("public object size exceeds AES-CTR stream range");
    }
    const std::uint64_t bytes_per_object =
        static_cast<std::uint64_t>(elements_per_object) * bytes_per_element;
    if (bytes_per_object != 0 &&
        object_index > (std::numeric_limits<std::uint64_t>::max)() / bytes_per_object) {
        throw std::overflow_error("public object offset exceeds AES-CTR stream range");
    }
    return object_index * bytes_per_object;
}

// A UniformRandomBitGenerator backed by buffered CSPRNG output. It is used
// only for secret randomness, never to reconstruct public values.
class SecureRng {
public:
    using result_type = std::uint64_t;

    static constexpr result_type min() { return 0; }
    static constexpr result_type max() { return (std::numeric_limits<result_type>::max)(); }

    SecureRng() = default;
    SecureRng(const SecureRng&) = delete;
    SecureRng& operator=(const SecureRng&) = delete;

    ~SecureRng() { secure_erase(buffer_.data(), buffer_.size()); }

    result_type operator()()
    {
        if (offset_ + sizeof(result_type) > buffer_.size()) {
            secure_random_bytes(buffer_.data(), buffer_.size());
            offset_ = 0;
        }

        result_type value = 0;
        std::memcpy(&value, buffer_.data() + offset_, sizeof(value));
        offset_ += sizeof(value);
        return value;
    }

private:
    std::array<std::uint8_t, 4096> buffer_{};
    std::size_t offset_ = buffer_.size();
};

} // namespace autohog
