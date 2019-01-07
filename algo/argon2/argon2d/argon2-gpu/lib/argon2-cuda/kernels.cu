/* For IDE: */
#ifndef __CUDACC__
#define __CUDACC__
#endif
#include <thread>
#include <chrono>


#include "argon2-cuda/kernels.h"
#include "argon2-gpu-common/argon2-common.h"

#include "argon2-cuda/cudaexception.h"

#include <stdexcept>
#ifndef NDEBUG
#include <iostream>

#endif

#define ARGON2_D  0
#define ARGON2_I  1
#define ARGON2_ID 2

#define ARGON2_VERSION_10 0x10
#define ARGON2_VERSION_13 0x13

#define ARGON2_BLOCK_SIZE 1024
#define ARGON2_DWORDS_IN_BLOCK (ARGON2_BLOCK_SIZE / 4)
#define ARGON2_QWORDS_IN_BLOCK (ARGON2_BLOCK_SIZE / 8)
#define ARGON2_SYNC_POINTS 4

#define THREADS_PER_LANE 32
#define QWORDS_PER_THREAD (ARGON2_QWORDS_IN_BLOCK / 32)

#define BLOCK_BYTES	32
#define OUT_BYTES	16
#define ARGON2_PREHASH_DIGEST_LENGTH	16
#define ARGON2_PREHASH_SEED_LENGTH		18

#define G(m, r, i, a, b, c, d) \
do { \
	a = a + b + m[blake2b_sigma[r][2 * i + 0]]; \
	d = rotr64(d ^ a, 32); \
	c = c + d; \
	b = rotr64(b ^ c, 24); \
	a = a + b + m[blake2b_sigma[r][2 * i + 1]]; \
	d = rotr64(d ^ a, 16); \
	c = c + d; \
	b = rotr64(b ^ c, 63); \
} while ((void)0, 0)

#define ROUND(m, v, r) \
do { \
	G(m, r, 0, v[0], v[4], v[ 8], v[12]); \
	G(m, r, 1, v[1], v[5], v[ 9], v[13]); \
	G(m, r, 2, v[2], v[6], v[10], v[14]); \
	G(m, r, 3, v[3], v[7], v[11], v[15]); \
	G(m, r, 4, v[0], v[5], v[10], v[15]); \
	G(m, r, 5, v[1], v[6], v[11], v[12]); \
	G(m, r, 6, v[2], v[7], v[ 8], v[13]); \
	G(m, r, 7, v[3], v[4], v[ 9], v[14]); \
} while ((void)0, 0)

namespace argon2 {
namespace cuda {

typedef struct blake2b_state_ {
	uint64_t h[8];
	uint64_t t[2];
	uint32_t buf[BLOCK_BYTES];
	uint32_t bufLen;
} blake2b_state;

__device__ uint64_t rotr64(uint64_t x, uint32_t n)
{
	return (x >> n) | (x << (64 - n));
}

__device__ void blake2b_init(blake2b_state *state, uint32_t outlen)
{
	uint64_t blake2b_IV[8] = {
			0x6A09E667F3BCC908, 0xBB67AE8584CAA73B,
			0x3C6EF372FE94F82B, 0xA54FF53A5F1D36F1,
			0x510E527FADE682D1, 0x9B05688C2B3E6C1F,
			0x1F83D9ABFB41BD6B, 0x5BE0CD19137E2179
	};

	state->t[1] = state->t[0] = 0;
	state->bufLen = 0;

	for(int i=0;i<8;i++) {
		state->h[i] = blake2b_IV[i];
	}

	state->h[0] ^= ((outlen * 4) | (1 << 16) | (1 << 24));
}

__device__ void blake2b_compress(blake2b_state *state, uint64_t *m, uint64_t f0)
{
	uint64_t v[16];

	uint64_t blake2b_IV[8] = {
			0x6A09E667F3BCC908, 0xBB67AE8584CAA73B,
			0x3C6EF372FE94F82B, 0xA54FF53A5F1D36F1,
			0x510E527FADE682D1, 0x9B05688C2B3E6C1F,
			0x1F83D9ABFB41BD6B, 0x5BE0CD19137E2179
	};
	unsigned int blake2b_sigma[12][16] = {
			{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15},
			{14, 10, 4, 8, 9, 15, 13, 6, 1, 12, 0, 2, 11, 7, 5, 3},
			{11, 8, 12, 0, 5, 2, 15, 13, 10, 14, 3, 6, 7, 1, 9, 4},
			{7, 9, 3, 1, 13, 12, 11, 14, 2, 6, 5, 10, 4, 0, 15, 8},
			{9, 0, 5, 7, 2, 4, 10, 15, 14, 1, 11, 12, 6, 8, 3, 13},
			{2, 12, 6, 10, 0, 11, 8, 3, 4, 13, 7, 5, 15, 14, 1, 9},
			{12, 5, 1, 15, 14, 13, 4, 10, 0, 7, 6, 3, 9, 2, 8, 11},
			{13, 11, 7, 14, 12, 1, 3, 9, 5, 0, 15, 4, 8, 6, 2, 10},
			{6, 15, 14, 9, 11, 3, 0, 8, 12, 2, 13, 7, 1, 4, 10, 5},
			{10, 2, 8, 4, 7, 6, 1, 5, 15, 11, 9, 14, 3, 12, 13, 0},
			{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15},
			{14, 10, 4, 8, 9, 15, 13, 6, 1, 12, 0, 2, 11, 7, 5, 3},
	};

	v[ 0] = state->h[0];
	v[ 1] = state->h[1];
	v[ 2] = state->h[2];
	v[ 3] = state->h[3];
	v[ 4] = state->h[4];
	v[ 5] = state->h[5];
	v[ 6] = state->h[6];
	v[ 7] = state->h[7];
	v[ 8] = blake2b_IV[0];
	v[ 9] = blake2b_IV[1];
	v[10] = blake2b_IV[2];
	v[11] = blake2b_IV[3];
	v[12] = blake2b_IV[4] ^ state->t[0];
	v[13] = blake2b_IV[5] ^ state->t[1];
	v[14] = blake2b_IV[6] ^ f0;
	v[15] = blake2b_IV[7];

	ROUND(m, v, 0);
	ROUND(m, v, 1);
	ROUND(m, v, 2);
	ROUND(m, v, 3);
	ROUND(m, v, 4);
	ROUND(m, v, 5);
	ROUND(m, v, 6);
	ROUND(m, v, 7);
	ROUND(m, v, 8);
	ROUND(m, v, 9);
	ROUND(m, v, 10);
	ROUND(m, v, 11);

	state->h[0] ^= v[0] ^ v[ 8];
	state->h[1] ^= v[1] ^ v[ 9];
	state->h[2] ^= v[2] ^ v[10];
	state->h[3] ^= v[3] ^ v[11];
	state->h[4] ^= v[4] ^ v[12];
	state->h[5] ^= v[5] ^ v[13];
	state->h[6] ^= v[6] ^ v[14];
	state->h[7] ^= v[7] ^ v[15];
}

__device__ void blake2b_incrementCounter(blake2b_state *state, uint32_t inc)
{
	state->t[0] += (inc * 4);
	state->t[1] += (state->t[0] < (inc * 4));
}

__device__ void blake2b_update(blake2b_state *state, uint32_t *in, int32_t inLen)
{
	if (state->bufLen + inLen > BLOCK_BYTES) {
		uint32_t temp[BLOCK_BYTES];
		uint32_t have = state->bufLen;
		uint32_t left = BLOCK_BYTES - have;

		for(int i=0;i<left;i++) {
			*(state->buf +  have + i) = in[i];
		}

		blake2b_incrementCounter(state, BLOCK_BYTES);
		blake2b_compress(state, (uint64_t*)state->buf, 0);

		state->bufLen = 0;
		inLen -= left;
		in += left;

		while (inLen > BLOCK_BYTES) {
			blake2b_incrementCounter(state, BLOCK_BYTES);

			for(int i=0;i<BLOCK_BYTES;i++) {
				temp[i] = in[i];
			}

			blake2b_compress(state, (uint64_t *)temp, 0);
			inLen -= BLOCK_BYTES;
			in += BLOCK_BYTES;
		}
	}
	for(int i=0;i<inLen;i++) {
		*(state->buf +  state->bufLen + i) = *(in + i);
	}
	state->bufLen += inLen;
}

__device__ void blake2b_final(blake2b_state *state, uint32_t *out, uint32_t outLen)
{
	blake2b_incrementCounter(state, state->bufLen);
	for(int i=0;i<BLOCK_BYTES - state->bufLen;i++) {
		*(state->buf + state->bufLen + i) = 0;
	}
	blake2b_compress(state, (uint64_t*)state->buf, 0xFFFFFFFFFFFFFFFF);
	for(int i=0;i<outLen;i++) {
		*(out + i) = *((uint32_t*)state->h + i);
	}
}

__device__ void blake2b_digestLong(uint32_t *out, uint32_t outLen,
								   uint32_t *in, uint32_t inLen)
{
	blake2b_state blake;

	if (outLen <= OUT_BYTES) {
		blake2b_init(&blake, outLen);

		blake.buf[0] = (outLen * 4);
		blake.bufLen = 1;

		blake2b_update(&blake, in, inLen);
		blake2b_final(&blake, out, outLen);
	} else {
		uint32_t out_buffer[OUT_BYTES];

		blake2b_init(&blake, OUT_BYTES);

		blake.buf[0] = (outLen * 4);
		blake.bufLen = 1;

		blake2b_update(&blake, in, inLen);
		blake2b_final(&blake, out_buffer, OUT_BYTES);

		for(int i=0;i<OUT_BYTES / 2;i++) {
			*(out + i) = *(out_buffer + i);
		}
		out += OUT_BYTES / 2;

		uint32_t toProduce = outLen - OUT_BYTES / 2;
		while (toProduce > OUT_BYTES) {
			blake2b_init(&blake, OUT_BYTES);
			blake2b_update(&blake, out_buffer, OUT_BYTES);
			blake2b_final(&blake, out_buffer, OUT_BYTES);

			for(int i=0;i<OUT_BYTES / 2;i++) {
				*(out + i) = *(out_buffer + i);
			}
			out += OUT_BYTES / 2;
			toProduce -= OUT_BYTES / 2;
		}

		blake2b_init(&blake, toProduce);
		blake2b_update(&blake, out_buffer, OUT_BYTES);
		blake2b_final(&blake, out, toProduce);
	}
}

__device__ uint64_t u64_build(uint32_t hi, uint32_t lo)
{
    return ((uint64_t)hi << 32) | (uint64_t)lo;
}

__device__ uint32_t u64_lo(uint64_t x)
{
    return (uint32_t)x;
}

__device__ uint32_t u64_hi(uint64_t x)
{
    return (uint32_t)(x >> 32);
}

struct u64_shuffle_buf {
    uint32_t lo[THREADS_PER_LANE];
    uint32_t hi[THREADS_PER_LANE];
};

__device__ uint64_t u64_shuffle(uint64_t v, uint32_t thread_src,
                                uint32_t thread, struct u64_shuffle_buf *buf)
{
    uint32_t lo = u64_lo(v);
    uint32_t hi = u64_hi(v);

    buf->lo[thread] = lo;
    buf->hi[thread] = hi;

    __syncthreads();

    lo = buf->lo[thread_src];
    hi = buf->hi[thread_src];

    return u64_build(hi, lo);
}

struct block_g {
    uint64_t data[ARGON2_QWORDS_IN_BLOCK];
};

struct block_th {
    uint64_t a, b, c, d;
};

__device__ uint64_t cmpeq_mask(uint32_t test, uint32_t ref)
{
    uint32_t x = -(uint32_t)(test == ref);
    return u64_build(x, x);
}

__device__ uint64_t block_th_get(const struct block_th *b, uint32_t idx)
{
    uint64_t res = 0;
    res ^= cmpeq_mask(idx, 0) & b->a;
    res ^= cmpeq_mask(idx, 1) & b->b;
    res ^= cmpeq_mask(idx, 2) & b->c;
    res ^= cmpeq_mask(idx, 3) & b->d;
    return res;
}

__device__ void block_th_set(struct block_th *b, uint32_t idx, uint64_t v)
{
    b->a ^= cmpeq_mask(idx, 0) & (v ^ b->a);
    b->b ^= cmpeq_mask(idx, 1) & (v ^ b->b);
    b->c ^= cmpeq_mask(idx, 2) & (v ^ b->c);
    b->d ^= cmpeq_mask(idx, 3) & (v ^ b->d);
}

__device__ void move_block(struct block_th *dst, const struct block_th *src)
{
    *dst = *src;
}

__device__ void xor_block(struct block_th *dst, const struct block_th *src)
{
    dst->a ^= src->a;
    dst->b ^= src->b;
    dst->c ^= src->c;
    dst->d ^= src->d;
}

__device__ void load_block(struct block_th *dst, const struct block_g *src,
                           uint32_t thread)
{
    dst->a = src->data[0 * THREADS_PER_LANE + thread];
    dst->b = src->data[1 * THREADS_PER_LANE + thread];
    dst->c = src->data[2 * THREADS_PER_LANE + thread];
    dst->d = src->data[3 * THREADS_PER_LANE + thread];
}

__device__ void load_block_xor(struct block_th *dst, const struct block_g *src,
                               uint32_t thread)
{
    dst->a ^= src->data[0 * THREADS_PER_LANE + thread];
    dst->b ^= src->data[1 * THREADS_PER_LANE + thread];
    dst->c ^= src->data[2 * THREADS_PER_LANE + thread];
    dst->d ^= src->data[3 * THREADS_PER_LANE + thread];
}

__device__ void store_block(struct block_g *dst, const struct block_th *src,
                            uint32_t thread)
{
    dst->data[0 * THREADS_PER_LANE + thread] = src->a;
    dst->data[1 * THREADS_PER_LANE + thread] = src->b;
    dst->data[2 * THREADS_PER_LANE + thread] = src->c;
    dst->data[3 * THREADS_PER_LANE + thread] = src->d;
}

__device__ uint64_t f(uint64_t x, uint64_t y)
{
    uint32_t xlo = u64_lo(x);
    uint32_t ylo = u64_lo(y);
    return x + y + 2 * u64_build(__umulhi(xlo, ylo), xlo * ylo);
}

__device__ void g(struct block_th *block)
{

	asm ("{"
	     ".reg .u64 d1, d2, a, b, c, d;\n\t"
		 ".reg .u32 s1, s2, s3, s4;\n\t"
		 "add.u64 d1, %0, %1;\n\t"
		 "cvt.u32.u64 s1, %0;\n\t"
		 "cvt.u32.u64 s2, %1;\n\t"
		 "mul.lo.u32 s3, s1, s2;\n\t"
		 "mul.hi.u32 s4, s1, s2;\n\t"
		 "mov.b64 a, {s3, s4};\n\t"
		 "shl.b64 d2, a, 1;\n\t"
		 "add.u64 a, d1, d2;\n\t"
		 "xor.b64 d1, %3, a;\n\t"
		 "mov.b64 {s1, s2}, d1;\n\t"
		 "mov.b64 d, {s2, s1};\n\t"
		 "add.u64 d1, %2, d;\n\t"
		 "cvt.u32.u64 s1, %2;\n\t"
		 "mul.lo.u32 s3, s1, s2;\n\t"
		 "mul.hi.u32 s4, s1, s2;\n\t"
		 "mov.b64 c, {s3, s4};\n\t"
		 "shl.b64 d2, c, 1;\n\t"
		 "add.u64 c, d1, d2;\n\t"
		 "xor.b64 d1, %1, c;\n\t"
		 "mov.b64 {s3, s4}, d1;\n\t"
		 "prmt.b32 s2, s3, s4, 0x6543;\n\t"
		 "prmt.b32 s1, s3, s4, 0x2107;\n\t"
		 "mov.b64 b, {s2, s1};\n\t"
		 "add.u64 d1, a, b;\n\t"
		 "cvt.u32.u64 s1, a;\n\t"
		 "mul.lo.u32 s3, s1, s2;\n\t"
		 "mul.hi.u32 s4, s1, s2;\n\t"
		 "mov.b64 a, {s3, s4};\n\t"
		 "shl.b64 d2, a, 1;\n\t"
		 "add.u64 %0, d1, d2;\n\t"
		 "xor.b64 d1, d, %0;\n\t"
		 "mov.b64 {s3, s4}, d1;\n\t"
		 "prmt.b32 s2, s3, s4, 0x5432;\n\t"
		 "prmt.b32 s1, s3, s4, 0x1076;\n\t"
		 "mov.b64 %3, {s2, s1};\n\t"
		 "add.u64 d1, c, %3;\n\t"
		 "cvt.u32.u64 s1, c;\n\t"
		 "mul.lo.u32 s3, s1, s2;\n\t"
		 "mul.hi.u32 s4, s1, s2;\n\t"
		 "mov.b64 c, {s3, s4};\n\t"
		 "shl.b64 d2, c, 1;\n\t"
		 "add.u64 %2, d1, d2;\n\t"
		 "xor.b64 d1, b, %2;\n\t"
		 "shl.b64 a, d1, 1;\n\t"
		 "shr.b64 b, d1, 63;\n\t"
		 "add.u64 %1, a, b;\n\t"
	"}" : "+l"(block->a), "+l"(block->b), "+l"(block->c), "+l"(block->d));

/*
    uint64_t a, b, c, d;
    a = block->a;
    b = block->b;
    c = block->c;
    d = block->d;

 	a = f(a, b);
    d = rotr64(d ^ a, 32);
    c = f(c, d);
    b = rotr64(b ^ c, 24);
    a = f(a, b);
    d = rotr64(d ^ a, 16);
    c = f(c, d);
    b = rotr64(b ^ c, 63);

    block->a = a;
    block->b = b;
    block->c = c;
    block->d = d; */
}

template<class shuffle>
__device__ void apply_shuffle(struct block_th *block, uint32_t thread,
                              struct u64_shuffle_buf *buf)
{
    for (uint32_t i = 0; i < QWORDS_PER_THREAD; i++) {
        uint32_t src_thr = shuffle::apply(thread, i);

        uint64_t v = block_th_get(block, i);
        v = u64_shuffle(v, src_thr, thread, buf);
        block_th_set(block, i, v);
    }
}

__device__ void transpose(struct block_th *block, uint32_t thread,
                          struct u64_shuffle_buf *buf)
{
    uint32_t thread_group = (thread & 0x0C) >> 2;
    for (uint32_t i = 1; i < QWORDS_PER_THREAD; i++) {
        uint32_t thr = (i << 2) ^ thread;
        uint32_t idx = thread_group ^ i;

        uint64_t v = block_th_get(block, idx);
        v = u64_shuffle(v, thr, thread, buf);
        block_th_set(block, idx, v);
    }
}

struct identity_shuffle {
    __device__ static uint32_t apply(uint32_t thread, uint32_t idx)
    {
        return thread;
    }
};

struct shift1_shuffle {
    __device__ static uint32_t apply(uint32_t thread, uint32_t idx)
    {
        return (thread & 0x1c) | ((thread + idx) & 0x3);
    }
};

struct unshift1_shuffle {
    __device__ static uint32_t apply(uint32_t thread, uint32_t idx)
    {
        idx = (QWORDS_PER_THREAD - idx) % QWORDS_PER_THREAD;

        return (thread & 0x1c) | ((thread + idx) & 0x3);
    }
};

struct shift2_shuffle {
    __device__ static uint32_t apply(uint32_t thread, uint32_t idx)
    {
        uint32_t lo = (thread & 0x1) | ((thread & 0x10) >> 3);
        lo = (lo + idx) & 0x3;
        return ((lo & 0x2) << 3) | (thread & 0xe) | (lo & 0x1);
    }
};

struct unshift2_shuffle {
    __device__ static uint32_t apply(uint32_t thread, uint32_t idx)
    {
        idx = (QWORDS_PER_THREAD - idx) % QWORDS_PER_THREAD;

        uint32_t lo = (thread & 0x1) | ((thread & 0x10) >> 3);
        lo = (lo + idx) & 0x3;
        return ((lo & 0x2) << 3) | (thread & 0xe) | (lo & 0x1);
    }
};

__device__ void shuffle_block(struct block_th *block, uint32_t thread,
                              struct u64_shuffle_buf *buf)
{
    transpose(block, thread, buf);

    g(block);

    apply_shuffle<shift1_shuffle>(block, thread, buf);

    g(block);

    apply_shuffle<unshift1_shuffle>(block, thread, buf);
    transpose(block, thread, buf);

    g(block);

    apply_shuffle<shift2_shuffle>(block, thread, buf);

    g(block);

    apply_shuffle<unshift2_shuffle>(block, thread, buf);
}

__device__ void next_addresses(struct block_th *addr, struct block_th *tmp,
                               uint32_t thread_input, uint32_t thread,
                               struct u64_shuffle_buf *buf)
{
    addr->a = u64_build(0, thread_input);
    addr->b = 0;
    addr->c = 0;
    addr->d = 0;

    shuffle_block(addr, thread, buf);

    addr->a ^= u64_build(0, thread_input);
    move_block(tmp, addr);

    shuffle_block(addr, thread, buf);

    xor_block(addr, tmp);
}

__device__ void compute_ref_pos(
        uint32_t lanes, uint32_t segment_blocks,
        uint32_t pass, uint32_t lane, uint32_t slice, uint32_t offset,
        uint32_t *ref_lane, uint32_t *ref_index)
{
    uint32_t lane_blocks = ARGON2_SYNC_POINTS * segment_blocks;

    *ref_lane = *ref_lane % lanes;

    uint32_t base;
    if (pass != 0) {
        base = lane_blocks - segment_blocks;
    } else {
        if (slice == 0) {
            *ref_lane = lane;
        }
        base = slice * segment_blocks;
    }

    uint32_t ref_area_size = base + offset - 1;
    if (*ref_lane != lane) {
        ref_area_size = min(ref_area_size, base);
    }

    *ref_index = __umulhi(*ref_index, *ref_index);
    *ref_index = ref_area_size - 1 - __umulhi(ref_area_size, *ref_index);

    if (pass != 0 && slice != ARGON2_SYNC_POINTS - 1) {
        *ref_index += (slice + 1) * segment_blocks;
        if (*ref_index >= lane_blocks) {
            *ref_index -= lane_blocks;
        }
    }
}

struct ref {
    uint32_t ref_lane;
    uint32_t ref_index;
};

/*
 * Refs hierarchy:
 * lanes -> passes -> slices -> blocks
 */
template<uint32_t type>
__global__ void argon2_precompute_kernel(
        struct ref *refs, uint32_t passes, uint32_t lanes,
        uint32_t segment_blocks)
{
    uint32_t block_id = blockIdx.y * blockDim.y + threadIdx.y;
    uint32_t warp = threadIdx.y;
    uint32_t thread = threadIdx.x;

    extern __shared__ struct u64_shuffle_buf shuffle_bufs[];
    struct u64_shuffle_buf *shuffle_buf = &shuffle_bufs[warp];

    uint32_t segment_addr_blocks = (segment_blocks + ARGON2_QWORDS_IN_BLOCK - 1)
            / ARGON2_QWORDS_IN_BLOCK;
    uint32_t block = block_id % segment_addr_blocks;
    uint32_t segment = block_id / segment_addr_blocks;

    uint32_t slice, pass, pass_id, lane;
    if (type == ARGON2_ID) {
        slice = segment % (ARGON2_SYNC_POINTS / 2);
        lane = segment / (ARGON2_SYNC_POINTS / 2);
        pass_id = pass = 0;
    } else {
        slice = segment % ARGON2_SYNC_POINTS;
        pass_id = segment / ARGON2_SYNC_POINTS;

        pass = pass_id % passes;
        lane = pass_id / passes;
    }

    struct block_th addr, tmp;

    uint32_t thread_input;
    switch (thread) {
    case 0:
        thread_input = pass;
        break;
    case 1:
        thread_input = lane;
        break;
    case 2:
        thread_input = slice;
        break;
    case 3:
        thread_input = lanes * segment_blocks * ARGON2_SYNC_POINTS;
        break;
    case 4:
        thread_input = passes;
        break;
    case 5:
        thread_input = type;
        break;
    case 6:
        thread_input = block + 1;
        break;
    default:
        thread_input = 0;
        break;
    }

    next_addresses(&addr, &tmp, thread_input, thread, shuffle_buf);

    refs += segment * segment_blocks;

    for (uint32_t i = 0; i < QWORDS_PER_THREAD; i++) {
        uint32_t pos = i * THREADS_PER_LANE + thread;
        uint32_t offset = block * ARGON2_QWORDS_IN_BLOCK + pos;
        if (offset < segment_blocks) {
            uint64_t v = block_th_get(&addr, i);
            uint32_t ref_index = u64_lo(v);
            uint32_t ref_lane  = u64_hi(v);

            compute_ref_pos(lanes, segment_blocks, pass, lane, slice, offset,
                            &ref_lane, &ref_index);

            refs[offset].ref_index = ref_index;
            refs[offset].ref_lane  = ref_lane;
        }
    }
}

template<uint32_t version>
__device__ void argon2_core(
        struct block_g *memory, struct block_g *mem_curr,
        struct block_th *prev, struct block_th *tmp,
        struct u64_shuffle_buf *shuffle_buf, uint32_t lanes,
        uint32_t thread, uint32_t pass, uint32_t ref_index, uint32_t ref_lane)
{
    struct block_g *mem_ref = memory + ref_index * lanes + ref_lane;

    if (version != ARGON2_VERSION_10 && pass != 0) {
        load_block(tmp, mem_curr, thread);
        load_block_xor(prev, mem_ref, thread);
        xor_block(tmp, prev);
    } else {
        load_block_xor(prev, mem_ref, thread);
        move_block(tmp, prev);
    }

    shuffle_block(prev, thread, shuffle_buf);

    xor_block(prev, tmp);

    store_block(mem_curr, prev, thread);
}

template<uint32_t type, uint32_t version>
__device__ void argon2_step_precompute(
        struct block_g *memory, struct block_g *mem_curr,
        struct block_th *prev, struct block_th *tmp,
        struct u64_shuffle_buf *shuffle_buf, const struct ref **refs,
        uint32_t lanes, uint32_t segment_blocks, uint32_t thread,
        uint32_t lane, uint32_t pass, uint32_t slice, uint32_t offset)
{
    uint32_t ref_index, ref_lane;
    if (type == ARGON2_I || (type == ARGON2_ID && pass == 0 &&
            slice < ARGON2_SYNC_POINTS / 2)) {
        ref_index = (*refs)->ref_index;
        ref_lane = (*refs)->ref_lane;
        (*refs)++;
    } else {
        uint64_t v = u64_shuffle(prev->a, 0, thread, shuffle_buf);
        ref_index = u64_lo(v);
        ref_lane  = u64_hi(v);

        compute_ref_pos(lanes, segment_blocks, pass, lane, slice, offset,
                        &ref_lane, &ref_index);
    }

    argon2_core<version>(memory, mem_curr, prev, tmp, shuffle_buf, lanes,
                         thread, pass, ref_index, ref_lane);
}

template<uint32_t type, uint32_t version>
__global__ void argon2_kernel_segment_precompute(
        struct block_g *memory, const struct ref *refs,
        uint32_t passes, uint32_t lanes, uint32_t segment_blocks,
        uint32_t pass, uint32_t slice)
{
    extern __shared__ struct u64_shuffle_buf shuffle_bufs[];
    struct u64_shuffle_buf *shuffle_buf =
            &shuffle_bufs[blockDim.y * threadIdx.z + threadIdx.y];

    uint32_t job_id = blockIdx.z * blockDim.z + threadIdx.z;
    uint32_t lane   = blockIdx.y * blockDim.y + threadIdx.y;
    uint32_t thread = threadIdx.x;

    uint32_t lane_blocks = ARGON2_SYNC_POINTS * segment_blocks;

    /* select job's memory region: */
    memory += (size_t)job_id * lanes * lane_blocks;

    struct block_th prev, tmp;

    struct block_g *mem_segment =
            memory + slice * segment_blocks * lanes + lane;
    struct block_g *mem_prev, *mem_curr;
    uint32_t start_offset = 0;
    if (pass == 0) {
        if (slice == 0) {
            mem_prev = mem_segment + 1 * lanes;
            mem_curr = mem_segment + 2 * lanes;
            start_offset = 2;
        } else {
            mem_prev = mem_segment - lanes;
            mem_curr = mem_segment;
        }
    } else {
        mem_prev = mem_segment + (slice == 0 ? lane_blocks * lanes : 0) - lanes;
        mem_curr = mem_segment;
    }

    load_block(&prev, mem_prev, thread);

    if (type == ARGON2_ID) {
        if (pass == 0 && slice < ARGON2_SYNC_POINTS / 2) {
            refs += lane * (lane_blocks / 2) + slice * segment_blocks;
            refs += start_offset;
        }
    } else {
        refs += (lane * passes + pass) * lane_blocks + slice * segment_blocks;
        refs += start_offset;
    }

    for (uint32_t offset = start_offset; offset < segment_blocks; ++offset) {
        argon2_step_precompute<type, version>(
                    memory, mem_curr, &prev, &tmp, shuffle_buf, &refs, lanes,
                    segment_blocks, thread, lane, pass, slice, offset);

        mem_curr += lanes;
    }
}

template<uint32_t type, uint32_t version>
__global__ void argon2_kernel_oneshot_precompute(
        struct block_g *memory, const struct ref *refs, uint32_t passes,
        uint32_t lanes, uint32_t segment_blocks)
{
    extern __shared__ struct u64_shuffle_buf shuffle_bufs[];
    struct u64_shuffle_buf *shuffle_buf =
            &shuffle_bufs[lanes * threadIdx.z + threadIdx.y];

    uint32_t job_id = blockIdx.z * blockDim.z + threadIdx.z;
    uint32_t lane   = threadIdx.y;
    uint32_t thread = threadIdx.x;

    uint32_t lane_blocks = ARGON2_SYNC_POINTS * segment_blocks;

    /* select job's memory region: */
    memory += (size_t)job_id * lanes * lane_blocks;

    struct block_th prev, tmp;

    struct block_g *mem_lane = memory + lane;
    struct block_g *mem_prev = mem_lane + 1 * lanes;
    struct block_g *mem_curr = mem_lane + 2 * lanes;

    load_block(&prev, mem_prev, thread);

    if (type == ARGON2_ID) {
        refs += lane * (lane_blocks / 2) + 2;
    } else {
        refs += lane * passes * lane_blocks + 2;
    }

    uint32_t skip = 2;
    for (uint32_t pass = 0; pass < passes; ++pass) {
        for (uint32_t slice = 0; slice < ARGON2_SYNC_POINTS; ++slice) {
            for (uint32_t offset = 0; offset < segment_blocks; ++offset) {
                if (skip > 0) {
                    --skip;
                    continue;
                }

                argon2_step_precompute<type, version>(
                            memory, mem_curr, &prev, &tmp, shuffle_buf, &refs,
                            lanes, segment_blocks, thread,
                            lane, pass, slice, offset);

                mem_curr += lanes;
            }

            __syncthreads();
        }

        mem_curr = mem_lane;
    }
}

template<uint32_t type, uint32_t version>
__device__ void argon2_step(
        struct block_g *memory, struct block_g *mem_curr,
        struct block_th *prev, struct block_th *tmp, struct block_th *addr,
        struct u64_shuffle_buf *shuffle_buf, uint32_t lanes,
        uint32_t segment_blocks, uint32_t thread, uint32_t *thread_input,
        uint32_t lane, uint32_t pass, uint32_t slice, uint32_t offset)
{
    uint32_t ref_index, ref_lane;

    if (type == ARGON2_I || (type == ARGON2_ID && pass == 0 &&
            slice < ARGON2_SYNC_POINTS / 2)) {
        uint32_t addr_index = offset % ARGON2_QWORDS_IN_BLOCK;
        if (addr_index == 0) {
            if (thread == 6) {
                ++*thread_input;
            }
            next_addresses(addr, tmp, *thread_input, thread, shuffle_buf);
        }

        uint32_t thr = addr_index % THREADS_PER_LANE;
        uint32_t idx = addr_index / THREADS_PER_LANE;

        uint64_t v = block_th_get(addr, idx);
        v = u64_shuffle(v, thr, thread, shuffle_buf);
        ref_index = u64_lo(v);
        ref_lane  = u64_hi(v);
    } else {
        uint64_t v = u64_shuffle(prev->a, 0, thread, shuffle_buf);
        ref_index = u64_lo(v);
        ref_lane  = u64_hi(v);
    }

    compute_ref_pos(lanes, segment_blocks, pass, lane, slice, offset,
                    &ref_lane, &ref_index);

    argon2_core<version>(memory, mem_curr, prev, tmp, shuffle_buf, lanes,
                         thread, pass, ref_index, ref_lane);
}

template<uint32_t type, uint32_t version>
__global__ void argon2_kernel_segment(
        struct block_g *memory, uint32_t passes, uint32_t lanes,
        uint32_t segment_blocks, uint32_t pass, uint32_t slice)
{
    extern __shared__ struct u64_shuffle_buf shuffle_bufs[];
    struct u64_shuffle_buf *shuffle_buf =
            &shuffle_bufs[blockDim.y * threadIdx.z + threadIdx.y];

    uint32_t job_id = blockIdx.z * blockDim.z + threadIdx.z;
    uint32_t lane   = blockIdx.y * blockDim.y + threadIdx.y;
    uint32_t thread = threadIdx.x;

    uint32_t lane_blocks = ARGON2_SYNC_POINTS * segment_blocks;

    /* select job's memory region: */
    memory += (size_t)job_id * lanes * lane_blocks;

    struct block_th prev, addr, tmp;
    uint32_t thread_input;

    if (type == ARGON2_I || type == ARGON2_ID) {
        switch (thread) {
        case 0:
            thread_input = pass;
            break;
        case 1:
            thread_input = lane;
            break;
        case 2:
            thread_input = slice;
            break;
        case 3:
            thread_input = lanes * lane_blocks;
            break;
        case 4:
            thread_input = passes;
            break;
        case 5:
            thread_input = type;
            break;
        default:
            thread_input = 0;
            break;
        }

        if (pass == 0 && slice == 0 && segment_blocks > 2) {
            if (thread == 6) {
                ++thread_input;
            }
            next_addresses(&addr, &tmp, thread_input, thread, shuffle_buf);
        }
    }

    struct block_g *mem_segment =
            memory + slice * segment_blocks * lanes + lane;
    struct block_g *mem_prev, *mem_curr;
    uint32_t start_offset = 0;
    if (pass == 0) {
        if (slice == 0) {
            mem_prev = mem_segment + 1 * lanes;
            mem_curr = mem_segment + 2 * lanes;
            start_offset = 2;
        } else {
            mem_prev = mem_segment - lanes;
            mem_curr = mem_segment;
        }
    } else {
        mem_prev = mem_segment + (slice == 0 ? lane_blocks * lanes : 0) - lanes;
        mem_curr = mem_segment;
    }

    load_block(&prev, mem_prev, thread);

    for (uint32_t offset = start_offset; offset < segment_blocks; ++offset) {
        argon2_step<type, version>(
                    memory, mem_curr, &prev, &tmp, &addr, shuffle_buf,
                    lanes, segment_blocks, thread, &thread_input,
                    lane, pass, slice, offset);

        mem_curr += lanes;
    }
}

template<uint32_t type, uint32_t version>
__global__ void argon2_kernel_oneshot(
        struct block_g *memory, uint32_t passes, uint32_t lanes,
        uint32_t segment_blocks)
{
    extern __shared__ struct u64_shuffle_buf shuffle_bufs[];
    struct u64_shuffle_buf *shuffle_buf =
            &shuffle_bufs[lanes * threadIdx.z + threadIdx.y];

    uint32_t job_id = blockIdx.z * blockDim.z + threadIdx.z;
    uint32_t lane   = threadIdx.y;
    uint32_t thread = threadIdx.x;

    uint32_t lane_blocks = ARGON2_SYNC_POINTS * segment_blocks;

    /* select job's memory region: */
    memory += (size_t)job_id * lanes * lane_blocks;

    struct block_th prev, addr, tmp;
    uint32_t thread_input;

    if (type == ARGON2_I || type == ARGON2_ID) {
        switch (thread) {
        case 1:
            thread_input = lane;
            break;
        case 3:
            thread_input = lanes * lane_blocks;
            break;
        case 4:
            thread_input = passes;
            break;
        case 5:
            thread_input = type;
            break;
        default:
            thread_input = 0;
            break;
        }

        if (segment_blocks > 2) {
            if (thread == 6) {
                ++thread_input;
            }
            next_addresses(&addr, &tmp, thread_input, thread, shuffle_buf);
        }
    }

    struct block_g *mem_lane = memory + lane;
    struct block_g *mem_prev = mem_lane + 1 * lanes;
    struct block_g *mem_curr = mem_lane + 2 * lanes;

    load_block(&prev, mem_prev, thread);

    uint32_t skip = 2;
    for (uint32_t pass = 0; pass < passes; ++pass) {
        for (uint32_t slice = 0; slice < ARGON2_SYNC_POINTS; ++slice) {
            for (uint32_t offset = 0; offset < segment_blocks; ++offset) {
                if (skip > 0) {
                    --skip;
                    continue;
                }

                argon2_step<type, version>(
                            memory, mem_curr, &prev, &tmp, &addr, shuffle_buf,
                            lanes, segment_blocks, thread, &thread_input,
                            lane, pass, slice, offset);

                mem_curr += lanes;
            }

            __syncthreads();

            if (type == ARGON2_I || type == ARGON2_ID) {
                if (thread == 2) {
                    ++thread_input;
                }
                if (thread == 6) {
                    thread_input = 0;
                }
            }
        }
        if (type == ARGON2_I) {
            if (thread == 0) {
                ++thread_input;
            }
            if (thread == 2) {
                thread_input = 0;
            }
        }
        mem_curr = mem_lane;
    }
}

__global__ void argon2_kernel_preseed(
		struct block_g *memory, uint32_t *seed, uint32_t lanes, uint32_t segment_blocks) {
	int job_id = blockIdx.x;
	int lane = threadIdx.x % lanes;
	int idx = threadIdx.x / lanes;

	/* select job's memory region: */
	memory += job_id * lanes * ARGON2_SYNC_POINTS * segment_blocks;
	seed += job_id * ARGON2_PREHASH_DIGEST_LENGTH;

	uint32_t initHash[ARGON2_PREHASH_SEED_LENGTH];
	for(int i=0;i<ARGON2_PREHASH_DIGEST_LENGTH;i++) {
		initHash[i] = seed[i];
	}

	initHash[ARGON2_PREHASH_DIGEST_LENGTH] = idx;
	initHash[ARGON2_PREHASH_DIGEST_LENGTH + 1] = lane;
	blake2b_digestLong((uint32_t*)(memory + lane + idx * lanes)->data, ARGON2_DWORDS_IN_BLOCK, initHash, ARGON2_PREHASH_SEED_LENGTH);
}

__global__ void argon2_kernel_finalize(
		struct block_g *memory, uint32_t *out, uint32_t outLen, uint32_t lanes, uint32_t segment_blocks) {
	int job_id = blockIdx.x;
	int thread = threadIdx.x;

	int lane_blocks = ARGON2_SYNC_POINTS * segment_blocks;
	/* select job's memory region: */
	memory += ((job_id + 1) * lanes * lane_blocks - lanes);
	out += job_id * outLen;
	struct block_g *dst = memory;

	for(int l=1;l<lanes;l++) {
		memory += 1;
		for (int i = 0; i < 4/*ARGON2_QWORDS_IN_BLOCK*/; i++) {
			dst->data[thread * 4 + i] ^= memory->data[thread * 4 + i];
		}
	}

	if(thread == 0) {
		blake2b_digestLong(out, outLen, (uint32_t *) dst, ARGON2_DWORDS_IN_BLOCK);
	}
}

KernelRunner::KernelRunner(uint32_t type, uint32_t version, uint32_t passes,
                           uint32_t lanes, uint32_t segmentBlocks,
                           size_t batchSize, size_t outLen, int32_t deviceIndex, bool bySegment, bool precompute)
    : type(type), version(version), passes(passes), lanes(lanes),
      segmentBlocks(segmentBlocks), batchSize(batchSize), outLen(outLen), bySegment(bySegment),
      precompute(precompute), stream(nullptr), memory(nullptr), seed(nullptr), seed_host(nullptr),
      out(nullptr), out_host(nullptr), refs(nullptr)
{

    int currentIndex = -1;
    CudaException::check(cudaGetDevice(&currentIndex));
    if (currentIndex != deviceIndex) {
        CudaException::check(cudaSetDevice(deviceIndex));
    }

    // FIXME: check overflow:
    size_t memorySize = batchSize * lanes * segmentBlocks
            * ARGON2_SYNC_POINTS * ARGON2_BLOCK_SIZE;

#ifndef NDEBUG
        std::cerr << "[INFO] Allocating " << memorySize << " bytes for memory..."
                  << std::endl;
#endif

    CudaException::check(cudaMalloc(&memory, memorySize));
    CudaException::check(cudaMalloc(&seed, batchSize * ARGON2_PREHASH_DIGEST_LENGTH * sizeof(uint32_t)));
    CudaException::check(cudaMalloc(&out, batchSize * outLen));
    CudaException::check(cudaMallocHost(&seed_host, batchSize * lanes * 2 * ARGON2_BLOCK_SIZE));
    CudaException::check(cudaMallocHost(&out_host, batchSize * outLen));

    CudaException::check(cudaStreamCreate(&stream));

    if ((type == ARGON2_I || type == ARGON2_ID) && precompute) {
        uint32_t segments =
                type == ARGON2_ID
                ? lanes * (ARGON2_SYNC_POINTS / 2)
                : passes * lanes * ARGON2_SYNC_POINTS;

        size_t refsSize = segments * segmentBlocks * sizeof(struct ref);

#ifndef NDEBUG
        std::cerr << "[INFO] Allocating " << refsSize << " bytes for refs..."
                  << std::endl;
#endif

        CudaException::check(cudaMalloc(&refs, refsSize));

        precomputeRefs();
        synchronize();
    }
}

void KernelRunner::precomputeRefs()
{
    struct ref *refs = (struct ref *)this->refs;

    uint32_t segmentAddrBlocks = (segmentBlocks + ARGON2_QWORDS_IN_BLOCK - 1)
            / ARGON2_QWORDS_IN_BLOCK;
    uint32_t segments =
            type == ARGON2_ID
            ? lanes * (ARGON2_SYNC_POINTS / 2)
            : passes * lanes * ARGON2_SYNC_POINTS;

    dim3 blocks = dim3(1, segments * segmentAddrBlocks);
    dim3 threads = dim3(THREADS_PER_LANE);

    size_t shmemSize = sizeof(struct u64_shuffle_buf);
    if (type == ARGON2_I) {
        argon2_precompute_kernel<ARGON2_I>
            <<<blocks, threads, shmemSize, stream>>>(
                refs, passes, lanes, segmentBlocks);
    } else {
        argon2_precompute_kernel<ARGON2_ID>
            <<<blocks, threads, shmemSize, stream>>>(
                refs, passes, lanes, segmentBlocks);
    }
}

KernelRunner::~KernelRunner()
{
    if (stream != nullptr) {
        cudaStreamDestroy(stream);
    }
    if (memory != nullptr) {
        cudaFree(memory);
    }
    if (refs != nullptr) {
        cudaFree(refs);
    }
    if (seed != nullptr) {
        cudaFree(seed);
    }
    if (out != nullptr) {
        cudaFree(out);
    }
    if (seed_host != nullptr) {
        cudaFreeHost(seed_host);
    }
    if (out_host != nullptr) {
        cudaFreeHost(out_host);
    }
}

void *KernelRunner::getSeedBuffer(int index) {
    return &((uint8_t*)seed_host)[index * ARGON2_PREHASH_DIGEST_LENGTH * sizeof(uint32_t)];
}

void *KernelRunner::getOutBuffer(int index) {
    return &((uint8_t*)out_host)[index * outLen];
}

void KernelRunner::writeInputMemory()
{
    CudaException::check(cudaMemcpyAsync(seed, seed_host, batchSize * ARGON2_PREHASH_DIGEST_LENGTH * sizeof(uint32_t),
                                         cudaMemcpyHostToDevice, stream));
}

void KernelRunner::readOutputMemory()
{
    CudaException::check(cudaMemcpyAsync(out_host, out, batchSize * outLen,
                                         cudaMemcpyDeviceToHost, stream));
}

void KernelRunner::runKernelSegment(uint32_t lanesPerBlock,
                                    uint32_t jobsPerBlock,
                                    uint32_t pass, uint32_t slice)
{
    if (lanesPerBlock > lanes || lanes % lanesPerBlock != 0) {
        throw std::logic_error("Invalid lanesPerBlock!");
    }

    if (jobsPerBlock > batchSize || batchSize % jobsPerBlock != 0) {
        throw std::logic_error("Invalid jobsPerBlock!");
    }

    struct block_g *memory_blocks = (struct block_g *)memory;
    dim3 blocks = dim3(1, lanes / lanesPerBlock, batchSize / jobsPerBlock);
    dim3 threads = dim3(THREADS_PER_LANE, lanesPerBlock, jobsPerBlock);
    uint32_t blockSize = lanesPerBlock * jobsPerBlock;
    uint32_t shared_size = blockSize * sizeof(struct u64_shuffle_buf);
    if (type == ARGON2_I) {
        if (precompute) {
            struct ref *refs = (struct ref *)this->refs;
            if (version == ARGON2_VERSION_10) {
                argon2_kernel_segment_precompute<ARGON2_I, ARGON2_VERSION_10>
                        <<<blocks, threads, shared_size, stream>>>(
                            memory_blocks, refs, passes, lanes, segmentBlocks,
                            pass, slice);
            } else {
                argon2_kernel_segment_precompute<ARGON2_I, ARGON2_VERSION_13>
                        <<<blocks, threads, shared_size, stream>>>(
                            memory_blocks, refs, passes, lanes, segmentBlocks,
                            pass, slice);
            }
        } else {
            if (version == ARGON2_VERSION_10) {
                argon2_kernel_segment<ARGON2_I, ARGON2_VERSION_10>
                        <<<blocks, threads, shared_size, stream>>>(
                            memory_blocks, passes, lanes, segmentBlocks,
                            pass, slice);
            } else {
                argon2_kernel_segment<ARGON2_I, ARGON2_VERSION_13>
                        <<<blocks, threads, shared_size, stream>>>(
                            memory_blocks, passes, lanes, segmentBlocks,
                            pass, slice);
            }
        }
    } else if (type == ARGON2_ID) {
        if (precompute) {
            struct ref *refs = (struct ref *)this->refs;
            if (version == ARGON2_VERSION_10) {
                argon2_kernel_segment_precompute<ARGON2_ID, ARGON2_VERSION_10>
                        <<<blocks, threads, shared_size, stream>>>(
                            memory_blocks, refs, passes, lanes, segmentBlocks,
                            pass, slice);
            } else {
                argon2_kernel_segment_precompute<ARGON2_ID, ARGON2_VERSION_13>
                        <<<blocks, threads, shared_size, stream>>>(
                            memory_blocks, refs, passes, lanes, segmentBlocks,
                            pass, slice);
            }
        } else {
            if (version == ARGON2_VERSION_10) {
                argon2_kernel_segment<ARGON2_ID, ARGON2_VERSION_10>
                        <<<blocks, threads, shared_size, stream>>>(
                            memory_blocks, passes, lanes, segmentBlocks,
                            pass, slice);
            } else {
                argon2_kernel_segment<ARGON2_ID, ARGON2_VERSION_13>
                        <<<blocks, threads, shared_size, stream>>>(
                            memory_blocks, passes, lanes, segmentBlocks,
                            pass, slice);
            }
        }
    } else {
        if (version == ARGON2_VERSION_10) {
            argon2_kernel_segment<ARGON2_D, ARGON2_VERSION_10>
                    <<<blocks, threads, shared_size, stream>>>(
                        memory_blocks, passes, lanes, segmentBlocks,
                        pass, slice);
        } else {
            argon2_kernel_segment<ARGON2_D, ARGON2_VERSION_13>
                    <<<blocks, threads, shared_size, stream>>>(
                        memory_blocks, passes, lanes, segmentBlocks,
                        pass, slice);
        }
    }
}

void KernelRunner::runKernelOneshot(uint32_t lanesPerBlock,
                                    uint32_t jobsPerBlock)
{
    if (lanesPerBlock != lanes) {
        throw std::logic_error("Invalid lanesPerBlock!");
    }

    if (jobsPerBlock > batchSize || batchSize % jobsPerBlock != 0) {
        throw std::logic_error("Invalid jobsPerBlock!");
    }

    struct block_g *memory_blocks = (struct block_g *)memory;
    dim3 blocks = dim3(1, 1, batchSize / jobsPerBlock);
    dim3 threads = dim3(THREADS_PER_LANE, lanes, jobsPerBlock);
    uint32_t blockSize = lanesPerBlock * jobsPerBlock;
    uint32_t shared_size = blockSize * sizeof(struct u64_shuffle_buf);
    if (type == ARGON2_I) {
        if (precompute) {
            struct ref *refs = (struct ref *)this->refs;
            if (version == ARGON2_VERSION_10) {
                argon2_kernel_oneshot_precompute<ARGON2_I, ARGON2_VERSION_10>
                        <<<blocks, threads, shared_size, stream>>>(
                            memory_blocks, refs, passes, lanes, segmentBlocks);
            } else {
                argon2_kernel_oneshot_precompute<ARGON2_I, ARGON2_VERSION_13>
                        <<<blocks, threads, shared_size, stream>>>(
                            memory_blocks, refs, passes, lanes, segmentBlocks);
            }
        } else {
            if (version == ARGON2_VERSION_10) {
                argon2_kernel_oneshot<ARGON2_I, ARGON2_VERSION_10>
                        <<<blocks, threads, shared_size, stream>>>(
                            memory_blocks, passes, lanes, segmentBlocks);
            } else {
                argon2_kernel_oneshot<ARGON2_I, ARGON2_VERSION_13>
                        <<<blocks, threads, shared_size, stream>>>(
                            memory_blocks, passes, lanes, segmentBlocks);
            }
        }
    } else if (type == ARGON2_ID) {
        if (precompute) {
            struct ref *refs = (struct ref *)this->refs;
            if (version == ARGON2_VERSION_10) {
                argon2_kernel_oneshot_precompute<ARGON2_ID, ARGON2_VERSION_10>
                        <<<blocks, threads, shared_size, stream>>>(
                            memory_blocks, refs, passes, lanes, segmentBlocks);
            } else {
                argon2_kernel_oneshot_precompute<ARGON2_ID, ARGON2_VERSION_13>
                        <<<blocks, threads, shared_size, stream>>>(
                            memory_blocks, refs, passes, lanes, segmentBlocks);
            }
        } else {
            if (version == ARGON2_VERSION_10) {
                argon2_kernel_oneshot<ARGON2_ID, ARGON2_VERSION_10>
                        <<<blocks, threads, shared_size, stream>>>(
                            memory_blocks, passes, lanes, segmentBlocks);
            } else {
                argon2_kernel_oneshot<ARGON2_ID, ARGON2_VERSION_13>
                        <<<blocks, threads, shared_size, stream>>>(
                            memory_blocks, passes, lanes, segmentBlocks);
            }
        }
    } else {
        if (version == ARGON2_VERSION_10) {
            argon2_kernel_oneshot<ARGON2_D, ARGON2_VERSION_10>
                    <<<blocks, threads, shared_size, stream>>>(
                        memory_blocks, passes, lanes, segmentBlocks);
        } else {
            argon2_kernel_oneshot<ARGON2_D, ARGON2_VERSION_13>
                    <<<blocks, threads, shared_size, stream>>>(
                        memory_blocks, passes, lanes, segmentBlocks);
        }
    }
}

void KernelRunner::runKernelPreseed() {
	struct block_g *memory_blocks = (struct block_g *)memory;
	argon2_kernel_preseed<<<batchSize, lanes * 2>>>(memory_blocks, (uint32_t *)seed, lanes, segmentBlocks);
}

void KernelRunner::runKernelFinalize() {
	struct block_g *memory_blocks = (struct block_g *)memory;
	argon2_kernel_finalize<<<batchSize, 32>>>(memory_blocks, (uint32_t *)out, outLen / 4, lanes, segmentBlocks);
}

void KernelRunner::run(uint32_t lanesPerBlock, uint32_t jobsPerBlock)
{
    timer = get_time();
    runKernelPreseed();
    if (bySegment) {
        for (uint32_t pass = 0; pass < passes; pass++) {
            for (uint32_t slice = 0; slice < ARGON2_SYNC_POINTS; slice++) {
                runKernelSegment(lanesPerBlock, jobsPerBlock, pass, slice);
            }
        }
    } else {
        runKernelOneshot(lanesPerBlock, jobsPerBlock);
    }
	runKernelFinalize();
}

uint64_t KernelRunner::finish()
{
	synchronize();
	return get_time() - timer;
}

void KernelRunner::synchronize() {
	while(cudaStreamQuery(stream) != cudaSuccess) {
		std::this_thread::sleep_for(std::chrono::milliseconds(1));
	};

	CudaException::check(cudaGetLastError());
}

uint64_t KernelRunner::get_time() {
    timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec * 1000000000 + t.tv_nsec;
}

} // cuda
} // argon2
