#include "argon2params.h"

#include "blake2b.h"

#include <cstring>
#include <algorithm>

#ifdef DEBUG
#include <cstdio>
#endif

namespace argon2 {

static void store32(void *dst, std::uint32_t v)
{
    auto out = static_cast<std::uint8_t *>(dst);
    *out++ = static_cast<std::uint8_t>(v); v >>= 8;
    *out++ = static_cast<std::uint8_t>(v); v >>= 8;
    *out++ = static_cast<std::uint8_t>(v); v >>= 8;
    *out++ = static_cast<std::uint8_t>(v);
}

Argon2Params::Argon2Params(
        std::size_t outLen,
        const void *salt, std::size_t saltLen,
        const void *secret, std::size_t secretLen,
        const void *ad, std::size_t adLen,
        std::size_t t_cost, std::size_t m_cost, std::size_t lanes)
    : salt(salt), secret(secret), ad(ad),
      outLen(outLen), saltLen(saltLen), secretLen(secretLen), adLen(adLen),
      t_cost(t_cost), m_cost(m_cost), lanes(lanes)
{
    // TODO validate inputs
    std::size_t segments = lanes * ARGON2_SYNC_POINTS;
    segmentBlocks = std::max(m_cost, 2 * segments) / segments;
}

void Argon2Params::digestLong(void *out, std::size_t outLen,
                              const void *in, std::size_t inLen)
{
    auto bout = static_cast<std::uint8_t *>(out);
    std::uint8_t outlen_bytes[sizeof(std::uint32_t)];
    Blake2b blake;

    store32(outlen_bytes, static_cast<std::uint32_t>(outLen));
    if (outLen <= Blake2b::OUT_BYTES) {
        blake.init(outLen);
        blake.update(outlen_bytes, sizeof(outlen_bytes));
        blake.update(in, inLen);
        blake.final(out, outLen);
    } else {
        std::uint8_t out_buffer[Blake2b::OUT_BYTES];

        blake.init(Blake2b::OUT_BYTES);
        blake.update(outlen_bytes, sizeof(outlen_bytes));
        blake.update(in, inLen);
        blake.final(out_buffer, Blake2b::OUT_BYTES);

        std::memcpy(bout, out_buffer, Blake2b::OUT_BYTES / 2);
        bout += Blake2b::OUT_BYTES / 2;

        std::size_t toProduce = outLen - Blake2b::OUT_BYTES / 2;
        while (toProduce > Blake2b::OUT_BYTES) {
            blake.init(Blake2b::OUT_BYTES);
            blake.update(out_buffer, Blake2b::OUT_BYTES);
            blake.final(out_buffer, Blake2b::OUT_BYTES);

            std::memcpy(bout, out_buffer, Blake2b::OUT_BYTES / 2);
            bout += Blake2b::OUT_BYTES / 2;
            toProduce -= Blake2b::OUT_BYTES / 2;
        }

        blake.init(toProduce);
        blake.update(out_buffer, Blake2b::OUT_BYTES);
        blake.final(bout, toProduce);
    }
}

void Argon2Params::initialHash(
        void *out, const void *pwd, std::size_t pwdLen,
        Type type, Version version) const
{
    Blake2b blake;
    std::uint8_t value[sizeof(std::uint32_t)];

    blake.init(ARGON2_PREHASH_DIGEST_LENGTH);

    store32(value, lanes);      blake.update(value, sizeof(value));
    store32(value, outLen);     blake.update(value, sizeof(value));
    store32(value, m_cost);     blake.update(value, sizeof(value));
    store32(value, t_cost);     blake.update(value, sizeof(value));
    store32(value, version);    blake.update(value, sizeof(value));
    store32(value, type);       blake.update(value, sizeof(value));
    store32(value, pwdLen);     blake.update(value, sizeof(value));
    blake.update(pwd, pwdLen);
    store32(value, saltLen);    blake.update(value, sizeof(value));
    blake.update(salt, saltLen);
    store32(value, secretLen);  blake.update(value, sizeof(value));
    blake.update(secret, secretLen);
    store32(value, adLen);      blake.update(value, sizeof(value));
    blake.update(ad, adLen);

    blake.final(out, ARGON2_PREHASH_DIGEST_LENGTH);
}

void Argon2Params::initialHashSameSalt(
        void *out, const void *pwd, std::size_t pwdLen,
        Type type, Version version) const
{
    Blake2b blake;
    std::uint8_t value[sizeof(std::uint32_t)];

    blake.init(ARGON2_PREHASH_DIGEST_LENGTH);

    store32(value, lanes);      blake.update(value, sizeof(value));
    store32(value, outLen);     blake.update(value, sizeof(value));
    store32(value, m_cost);     blake.update(value, sizeof(value));
    store32(value, t_cost);     blake.update(value, sizeof(value));
    store32(value, version);    blake.update(value, sizeof(value));
    store32(value, type);       blake.update(value, sizeof(value));
    store32(value, pwdLen);     blake.update(value, sizeof(value));
    blake.update(pwd, pwdLen);
    store32(value, pwdLen);    blake.update(value, sizeof(value));
    blake.update(pwd, pwdLen);
    store32(value, secretLen);  blake.update(value, sizeof(value));
    blake.update(secret, secretLen);
    store32(value, adLen);      blake.update(value, sizeof(value));
    blake.update(ad, adLen);

    blake.final(out, ARGON2_PREHASH_DIGEST_LENGTH);
}


void Argon2Params::fillFirstBlocks(
        void *memory, const void *pwd, std::size_t pwdLen,
        Type type, Version version) const
{
	auto bmemory = static_cast<std::uint8_t *>(memory);
	initialHash(bmemory, pwd, pwdLen, type, version);

#ifdef DEBUG
    std::fprintf(stderr, "Initial hash: ");
    for (std::size_t i = 0; i < ARGON2_PREHASH_DIGEST_LENGTH; i++) {
        std::fprintf(stderr, "%02x", (unsigned int)bmemory[i]);
    }
    std::fprintf(stderr, "\n");
#endif
}

void Argon2Params::fillFirstBlocksSameSalt(
        void *memory, const void *pwd, std::size_t pwdLen,
        Type type, Version version) const
{
	auto bmemory = static_cast<std::uint8_t *>(memory);
    initialHashSameSalt(bmemory, pwd, pwdLen, type, version);

#ifdef DEBUG
    std::fprintf(stderr, "Initial hash: ");
    for (std::size_t i = 0; i < ARGON2_PREHASH_DIGEST_LENGTH; i++) {
        std::fprintf(stderr, "%02x", (unsigned int)bmemory[i]);
    }
    std::fprintf(stderr, "\n");
#endif
}


void Argon2Params::finalize(void *out, const void *memory) const
{
    digestLong(out, outLen, memory, ARGON2_BLOCK_SIZE);
}

} // namespace argon2

