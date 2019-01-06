/* C compatibility For dumb IDEs: */
#ifndef __OPENCL_VERSION__
#ifndef __cplusplus
typedef int bool;
#endif
typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int uint;
typedef unsigned long ulong;
typedef unsigned long size_t;
typedef long ptrdiff_t;
typedef size_t uintptr_t;
typedef ptrdiff_t intptr_t;
#ifndef __kernel
#define __kernel
#endif
#ifndef __global
#define __global
#endif
#ifndef __private
#define __private
#endif
#ifndef __local
#define __local
#endif
#ifndef __constant
#define __constant const
#endif
#endif /* __OPENCL_VERSION__ */

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

#ifndef ARGON2_VERSION
#define ARGON2_VERSION ARGON2_VERSION_13
#endif

#ifndef ARGON2_TYPE
#define ARGON2_TYPE ARGON2_I
#endif

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

typedef struct blake2b_state_ {
    ulong h[8];
    ulong t[2];
    uint buf[BLOCK_BYTES];
    uint bufLen;
} blake2b_state;

ulong rotr64(ulong x, ulong n)
{
	return rotate(x, 64 - n);
}

void blake2b_init(blake2b_state *state, uint outlen)
{
    ulong blake2b_IV[8] = {
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

void blake2b_compress(blake2b_state *state, ulong *m, ulong f0)
{
    ulong v[16];

    ulong blake2b_IV[8] = {
            0x6A09E667F3BCC908, 0xBB67AE8584CAA73B,
            0x3C6EF372FE94F82B, 0xA54FF53A5F1D36F1,
            0x510E527FADE682D1, 0x9B05688C2B3E6C1F,
            0x1F83D9ABFB41BD6B, 0x5BE0CD19137E2179
    };
    uint blake2b_sigma[12][16] = {
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

void blake2b_incrementCounter(blake2b_state *state, uint inc)
{
    state->t[0] += (inc * 4);
    state->t[1] += (state->t[0] < (inc * 4));
}

void blake2b_update(blake2b_state *state, uint *in, int inLen)
{
    if (state->bufLen + inLen > BLOCK_BYTES) {
        uint temp[BLOCK_BYTES];
        uint have = state->bufLen;
        uint left = BLOCK_BYTES - have;

        for(int i=0;i<left;i++) {
            *(state->buf +  have + i) = in[i];
        }

        blake2b_incrementCounter(state, BLOCK_BYTES);
        blake2b_compress(state, (ulong*)state->buf, 0);

        state->bufLen = 0;
        inLen -= left;
        in += left;

        while (inLen > BLOCK_BYTES) {
            blake2b_incrementCounter(state, BLOCK_BYTES);

            for(int i=0;i<BLOCK_BYTES;i++) {
                temp[i] = in[i];
            }

            blake2b_compress(state, (ulong*)temp, 0);
            inLen -= BLOCK_BYTES;
            in += BLOCK_BYTES;
        }
    }
    for(int i=0;i<inLen;i++) {
        *(state->buf +  state->bufLen + i) = *(in + i);
    }
    state->bufLen += inLen;
}

void blake2b_final(blake2b_state *state, uint *out, uint outLen)
{
    blake2b_incrementCounter(state, state->bufLen);
    for(int i=0;i<BLOCK_BYTES - state->bufLen;i++) {
        *(state->buf + state->bufLen + i) = 0;
    }
    blake2b_compress(state, (ulong*)state->buf, 0xFFFFFFFFFFFFFFFF);
    for(int i=0;i<outLen;i++) {
        *(out + i) = *((uint*)state->h + i);
    }
}

void blake2b_update_global(blake2b_state *state, __global uint *in, int inLen)
{
	if (state->bufLen + inLen > BLOCK_BYTES) {
		uint temp[BLOCK_BYTES];
		uint have = state->bufLen;
		uint left = BLOCK_BYTES - have;

		for(int i=0;i<left;i++) {
			*(state->buf +  have + i) = in[i];
		}

		blake2b_incrementCounter(state, BLOCK_BYTES);
		blake2b_compress(state, (ulong*)state->buf, 0);

		state->bufLen = 0;
		inLen -= left;
		in += left;

		while (inLen > BLOCK_BYTES) {
			blake2b_incrementCounter(state, BLOCK_BYTES);

			for(int i=0;i<BLOCK_BYTES;i++) {
				temp[i] = in[i];
			}

			blake2b_compress(state, (ulong*)temp, 0);
			inLen -= BLOCK_BYTES;
			in += BLOCK_BYTES;
		}
	}
	for(int i=0;i<inLen;i++) {
		*(state->buf +  state->bufLen + i) = *(in + i);
	}
	state->bufLen += inLen;
}

void blake2b_final_global(blake2b_state *state, __global uint *out, uint outLen)
{
	blake2b_incrementCounter(state, state->bufLen);
	for(int i=0;i<BLOCK_BYTES - state->bufLen;i++) {
		*(state->buf + state->bufLen + i) = 0;
	}
	blake2b_compress(state, (ulong*)state->buf, 0xFFFFFFFFFFFFFFFF);
	for(int i=0;i<outLen;i++) {
		*(out + i) = *((uint*)state->h + i);
	}
}

void blake2b_digestLong(__global uint *out, uint outLen,
                                   __global uint *in, uint inLen)
{
    blake2b_state blake;

    if (outLen <= OUT_BYTES) {
        blake2b_init(&blake, outLen);

        blake.buf[0] = (outLen * 4);
        blake.bufLen = 1;

        blake2b_update_global(&blake, in, inLen);
        blake2b_final_global(&blake, out, outLen);
    } else {
        uint out_buffer[OUT_BYTES];

        blake2b_init(&blake, OUT_BYTES);

        blake.buf[0] = (outLen * 4);
        blake.bufLen = 1;

        blake2b_update_global(&blake, in, inLen);
        blake2b_final(&blake, out_buffer, OUT_BYTES);

        for(int i=0;i<OUT_BYTES / 2;i++) {
            *(out + i) = *(out_buffer + i);
        }
        out += OUT_BYTES / 2;

        uint toProduce = outLen - OUT_BYTES / 2;
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
        blake2b_final_global(&blake, out, toProduce);
    }
}

ulong u64_build(uint hi, uint lo)
{
    return upsample(hi, lo);
}

uint u64_lo(ulong x)
{
    return (uint)x;
}

uint u64_hi(ulong x)
{
    return (uint)(x >> 32);
}

struct u64_shuffle_buf {
    uint lo[THREADS_PER_LANE];
    uint hi[THREADS_PER_LANE];
};

ulong u64_shuffle(ulong v, uint thread_src, uint thread,
                  __local struct u64_shuffle_buf *buf)
{
    uint lo = u64_lo(v);
    uint hi = u64_hi(v);

    buf->lo[thread] = lo;
    buf->hi[thread] = hi;

    barrier(CLK_LOCAL_MEM_FENCE);

    lo = buf->lo[thread_src];
    hi = buf->hi[thread_src];

    return u64_build(hi, lo);
}

struct block_g {
    ulong data[ARGON2_QWORDS_IN_BLOCK];
};

struct block_th {
    ulong a, b, c, d;
};

ulong cmpeq_mask(uint test, uint ref)
{
    uint x = -(uint)(test == ref);
    return u64_build(x, x);
}

ulong block_th_get(const struct block_th *b, uint idx)
{
    ulong res = 0;
    res ^= cmpeq_mask(idx, 0) & b->a;
    res ^= cmpeq_mask(idx, 1) & b->b;
    res ^= cmpeq_mask(idx, 2) & b->c;
    res ^= cmpeq_mask(idx, 3) & b->d;
    return res;
}

void block_th_set(struct block_th *b, uint idx, ulong v)
{
    b->a ^= cmpeq_mask(idx, 0) & (v ^ b->a);
    b->b ^= cmpeq_mask(idx, 1) & (v ^ b->b);
    b->c ^= cmpeq_mask(idx, 2) & (v ^ b->c);
    b->d ^= cmpeq_mask(idx, 3) & (v ^ b->d);
}

void move_block(struct block_th *dst, const struct block_th *src)
{
    *dst = *src;
}

void xor_block(struct block_th *dst, const struct block_th *src)
{
    dst->a ^= src->a;
    dst->b ^= src->b;
    dst->c ^= src->c;
    dst->d ^= src->d;
}

void load_block(struct block_th *dst, __global const struct block_g *src,
                uint thread)
{
    dst->a = src->data[0 * THREADS_PER_LANE + thread];
    dst->b = src->data[1 * THREADS_PER_LANE + thread];
    dst->c = src->data[2 * THREADS_PER_LANE + thread];
    dst->d = src->data[3 * THREADS_PER_LANE + thread];
}

void load_block_xor(struct block_th *dst, __global const struct block_g *src,
                    uint thread)
{
    dst->a ^= src->data[0 * THREADS_PER_LANE + thread];
    dst->b ^= src->data[1 * THREADS_PER_LANE + thread];
    dst->c ^= src->data[2 * THREADS_PER_LANE + thread];
    dst->d ^= src->data[3 * THREADS_PER_LANE + thread];
}

void store_block(__global struct block_g *dst, const struct block_th *src,
                 uint thread)
{
    dst->data[0 * THREADS_PER_LANE + thread] = src->a;
    dst->data[1 * THREADS_PER_LANE + thread] = src->b;
    dst->data[2 * THREADS_PER_LANE + thread] = src->c;
    dst->data[3 * THREADS_PER_LANE + thread] = src->d;
}

ulong f(ulong x, ulong y)
{
    uint xlo = u64_lo(x);
    uint ylo = u64_lo(y);
    return x + y + 2 * u64_build(mul_hi(xlo, ylo), xlo * ylo);
}

void g(struct block_th *block)
{
    ulong a, b, c, d;
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
    block->d = d;
}

uint apply_shuffle_shift1(uint thread, uint idx)
{
    return (thread & 0x1c) | ((thread + idx) & 0x3);
}

uint apply_shuffle_unshift1(uint thread, uint idx)
{
    idx = (QWORDS_PER_THREAD - idx) % QWORDS_PER_THREAD;

    return apply_shuffle_shift1(thread, idx);
}

uint apply_shuffle_shift2(uint thread, uint idx)
{
    uint lo = (thread & 0x1) | ((thread & 0x10) >> 3);
    lo = (lo + idx) & 0x3;
    return ((lo & 0x2) << 3) | (thread & 0xe) | (lo & 0x1);
}

uint apply_shuffle_unshift2(uint thread, uint idx)
{
    idx = (QWORDS_PER_THREAD - idx) % QWORDS_PER_THREAD;

    return apply_shuffle_shift2(thread, idx);
}

void shuffle_shift1(struct block_th *block, uint thread,
                    __local struct u64_shuffle_buf *buf)
{
    for (uint i = 0; i < QWORDS_PER_THREAD; i++) {
        uint src_thr = apply_shuffle_shift1(thread, i);

        ulong v = block_th_get(block, i);
        v = u64_shuffle(v, src_thr, thread, buf);
        block_th_set(block, i, v);
    }
}

void shuffle_unshift1(struct block_th *block, uint thread,
                      __local struct u64_shuffle_buf *buf)
{
    for (uint i = 0; i < QWORDS_PER_THREAD; i++) {
        uint src_thr = apply_shuffle_unshift1(thread, i);

        ulong v = block_th_get(block, i);
        v = u64_shuffle(v, src_thr, thread, buf);
        block_th_set(block, i, v);
    }
}

void shuffle_shift2(struct block_th *block, uint thread,
                    __local struct u64_shuffle_buf *buf)
{
    for (uint i = 0; i < QWORDS_PER_THREAD; i++) {
        uint src_thr = apply_shuffle_shift2(thread, i);

        ulong v = block_th_get(block, i);
        v = u64_shuffle(v, src_thr, thread, buf);
        block_th_set(block, i, v);
    }
}

void shuffle_unshift2(struct block_th *block, uint thread,
                      __local struct u64_shuffle_buf *buf)
{
    for (uint i = 0; i < QWORDS_PER_THREAD; i++) {
        uint src_thr = apply_shuffle_unshift2(thread, i);

        ulong v = block_th_get(block, i);
        v = u64_shuffle(v, src_thr, thread, buf);
        block_th_set(block, i, v);
    }
}

void transpose(struct block_th *block, uint thread,
               __local struct u64_shuffle_buf *buf)
{
    uint thread_group = (thread & 0x0C) >> 2;
    for (uint i = 1; i < QWORDS_PER_THREAD; i++) {
        uint thr = (i << 2) ^ thread;
        uint idx = thread_group ^ i;

        ulong v = block_th_get(block, idx);
        v = u64_shuffle(v, thr, thread, buf);
        block_th_set(block, idx, v);
    }
}

void shuffle_block(struct block_th *block, uint thread,
                   __local struct u64_shuffle_buf *buf)
{
    transpose(block, thread, buf);

    g(block);

    shuffle_shift1(block, thread, buf);

    g(block);

    shuffle_unshift1(block, thread, buf);
    transpose(block, thread, buf);

    g(block);

    shuffle_shift2(block, thread, buf);

    g(block);

    shuffle_unshift2(block, thread, buf);
}

void compute_ref_pos(uint lanes, uint segment_blocks,
                     uint pass, uint lane, uint slice, uint offset,
                     uint *ref_lane, uint *ref_index)
{
    uint lane_blocks = ARGON2_SYNC_POINTS * segment_blocks;

    *ref_lane = *ref_lane % lanes;

    uint base;
    if (pass != 0) {
        base = lane_blocks - segment_blocks;
    } else {
        if (slice == 0) {
            *ref_lane = lane;
        }
        base = slice * segment_blocks;
    }

    uint ref_area_size = base + offset - 1;
    if (*ref_lane != lane) {
        ref_area_size = min(ref_area_size, base);
    }

    *ref_index = mul_hi(*ref_index, *ref_index);
    *ref_index = ref_area_size - 1 - mul_hi(ref_area_size, *ref_index);

    if (pass != 0 && slice != ARGON2_SYNC_POINTS - 1) {
        *ref_index += (slice + 1) * segment_blocks;
        if (*ref_index >= lane_blocks) {
            *ref_index -= lane_blocks;
        }
    }
}

void argon2_core(
        __global struct block_g *memory, __global struct block_g *mem_curr,
        struct block_th *prev, struct block_th *tmp,
        __local struct u64_shuffle_buf *shuffle_buf, uint lanes,
        uint thread, uint pass, uint ref_index, uint ref_lane)
{
    __global struct block_g *mem_ref;
    mem_ref = memory + ref_index * lanes + ref_lane;

#if ARGON2_VERSION == ARGON2_VERSION_10
    load_block_xor(prev, mem_ref, thread);
    move_block(tmp, prev);
#else
    if (pass != 0) {
        load_block(tmp, mem_curr, thread);
        load_block_xor(prev, mem_ref, thread);
        xor_block(tmp, prev);
    } else {
        load_block_xor(prev, mem_ref, thread);
        move_block(tmp, prev);
    }
#endif

    shuffle_block(prev, thread, shuffle_buf);

    xor_block(prev, tmp);

    store_block(mem_curr, prev, thread);
}

void next_addresses(struct block_th *addr, struct block_th *tmp,
                    uint thread_input, uint thread,
                    __local struct u64_shuffle_buf *buf)
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

#if ARGON2_TYPE == ARGON2_I || ARGON2_TYPE == ARGON2_ID
struct ref {
    uint ref_lane;
    uint ref_index;
};

/*
 * Refs hierarchy:
 * lanes -> passes -> slices -> blocks
 */
__kernel void argon2_precompute_kernel(
        __local struct u64_shuffle_buf *shuffle_bufs, __global struct ref *refs,
        uint passes, uint lanes, uint segment_blocks)
{
    uint block_id = get_global_id(0) / THREADS_PER_LANE;
    uint warp = get_local_id(0) / THREADS_PER_LANE;
    uint thread = get_local_id(0) % THREADS_PER_LANE;

    __local struct u64_shuffle_buf *shuffle_buf = &shuffle_bufs[warp];

    uint segment_addr_blocks = (segment_blocks + ARGON2_QWORDS_IN_BLOCK - 1)
            / ARGON2_QWORDS_IN_BLOCK;
    uint block = block_id % segment_addr_blocks;
    uint segment = block_id / segment_addr_blocks;

    uint slice, pass, lane;
#if ARGON2_TYPE == ARGON2_ID
    slice = segment % (ARGON2_SYNC_POINTS / 2);
    lane = segment / (ARGON2_SYNC_POINTS / 2);
    pass = 0;
#else
    uint pass_id;

    slice = segment % ARGON2_SYNC_POINTS;
    pass_id = segment / ARGON2_SYNC_POINTS;

    pass = pass_id % passes;
    lane = pass_id / passes;
#endif

    struct block_th addr, tmp;

    uint thread_input;
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
        thread_input = ARGON2_TYPE;
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

    for (uint i = 0; i < QWORDS_PER_THREAD; i++) {
        uint pos = i * THREADS_PER_LANE + thread;
        uint offset = block * ARGON2_QWORDS_IN_BLOCK + pos;
        if (offset < segment_blocks) {
            ulong v = block_th_get(&addr, i);
            uint ref_index = u64_lo(v);
            uint ref_lane  = u64_hi(v);

            compute_ref_pos(lanes, segment_blocks, pass, lane, slice, offset,
                            &ref_lane, &ref_index);

            refs[offset].ref_index = ref_index;
            refs[offset].ref_lane  = ref_lane;
        }
    }
}

void argon2_step_precompute(
        __global struct block_g *memory, __global struct block_g *mem_curr,
        struct block_th *prev, struct block_th *tmp,
        __local struct u64_shuffle_buf *shuffle_buf,
        __global const struct ref **refs,
        uint lanes, uint segment_blocks, uint thread,
        uint lane, uint pass, uint slice, uint offset)
{
    uint ref_index, ref_lane;
    bool data_independent;
#if ARGON2_TYPE == ARGON2_I
    data_independent = true;
#elif ARGON2_TYPE == ARGON2_ID
    data_independent = pass == 0 && slice < ARGON2_SYNC_POINTS / 2;
#else
    data_independent = false;
#endif
    if (data_independent) {
        ref_index = (*refs)->ref_index;
        ref_lane = (*refs)->ref_lane;
        (*refs)++;
    } else {
        ulong v = u64_shuffle(prev->a, 0, thread, shuffle_buf);
        ref_index = u64_lo(v);
        ref_lane  = u64_hi(v);

        compute_ref_pos(lanes, segment_blocks, pass, lane, slice, offset,
                        &ref_lane, &ref_index);
    }

    argon2_core(memory, mem_curr, prev, tmp, shuffle_buf, lanes, thread, pass,
                ref_index, ref_lane);
}

__kernel void argon2_kernel_segment_precompute(
        __local struct u64_shuffle_buf *shuffle_bufs,
        __global struct block_g *memory, __global const struct ref *refs,
        uint passes, uint lanes, uint segment_blocks,
        uint pass, uint slice)
{
    uint job_id = get_global_id(1);
    uint lane   = get_global_id(0) / THREADS_PER_LANE;
    uint warp   = (get_local_id(1) * get_local_size(0) + get_local_id(0))
            / THREADS_PER_LANE;
    uint thread = get_local_id(0) % THREADS_PER_LANE;

    __local struct u64_shuffle_buf *shuffle_buf = &shuffle_bufs[warp];

    uint lane_blocks = ARGON2_SYNC_POINTS * segment_blocks;

    /* select job's memory region: */
    memory += (size_t)job_id * lanes * lane_blocks;

    struct block_th prev, tmp;

    __global struct block_g *mem_segment =
            memory + slice * segment_blocks * lanes + lane;
    __global struct block_g *mem_prev, *mem_curr;
    uint start_offset = 0;
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

#if ARGON2_TYPE == ARGON2_ID
        if (pass == 0 && slice < ARGON2_SYNC_POINTS / 2) {
            refs += lane * (lane_blocks / 2) + slice * segment_blocks;
            refs += start_offset;
        }
#else
        refs += (lane * passes + pass) * lane_blocks + slice * segment_blocks;
        refs += start_offset;
#endif

    for (uint offset = start_offset; offset < segment_blocks; ++offset) {
        argon2_step_precompute(
                    memory, mem_curr, &prev, &tmp, shuffle_buf, &refs, lanes,
                    segment_blocks, thread, lane, pass, slice, offset);

        mem_curr += lanes;
    }
}

__kernel void argon2_kernel_oneshot_precompute(
        __local struct u64_shuffle_buf *shuffle_bufs,
        __global struct block_g *memory, __global const struct ref *refs,
        uint passes, uint lanes, uint segment_blocks)
{
    uint job_id = get_global_id(1);
    uint lane   = get_global_id(0) / THREADS_PER_LANE;
    uint warp   = get_local_id(1) * lanes + get_local_id(0) / THREADS_PER_LANE;
    uint thread = get_local_id(0) % THREADS_PER_LANE;

    __local struct u64_shuffle_buf *shuffle_buf = &shuffle_bufs[warp];

    uint lane_blocks = ARGON2_SYNC_POINTS * segment_blocks;

    /* select job's memory region: */
    memory += (size_t)job_id * lanes * lane_blocks;

    struct block_th prev, tmp;

    __global struct block_g *mem_lane = memory + lane;
    __global struct block_g *mem_prev = mem_lane + 1 * lanes;
    __global struct block_g *mem_curr = mem_lane + 2 * lanes;

    load_block(&prev, mem_prev, thread);

#if ARGON2_TYPE == ARGON2_ID
    refs += lane * (lane_blocks / 2) + 2;
#else
    refs += lane * passes * lane_blocks + 2;
#endif

    uint skip = 2;
    for (uint pass = 0; pass < passes; ++pass) {
        for (uint slice = 0; slice < ARGON2_SYNC_POINTS; ++slice) {
            for (uint offset = 0; offset < segment_blocks; ++offset) {
                if (skip > 0) {
                    --skip;
                    continue;
                }

                argon2_step_precompute(
                            memory, mem_curr, &prev, &tmp, shuffle_buf, &refs,
                            lanes, segment_blocks, thread,
                            lane, pass, slice, offset);

                mem_curr += lanes;
            }

            barrier(CLK_GLOBAL_MEM_FENCE);
        }

        mem_curr = mem_lane;
    }
}
#endif /* ARGON2_TYPE == ARGON2_I || ARGON2_TYPE == ARGON2_ID */

void argon2_step(
        __global struct block_g *memory, __global struct block_g *mem_curr,
        struct block_th *prev, struct block_th *tmp, struct block_th *addr,
        __local struct u64_shuffle_buf *shuffle_buf,
        uint lanes, uint segment_blocks, uint thread, uint *thread_input,
        uint lane, uint pass, uint slice, uint offset)
{
    uint ref_index, ref_lane;
    bool data_independent;
#if ARGON2_TYPE == ARGON2_I
    data_independent = true;
#elif ARGON2_TYPE == ARGON2_ID
    data_independent = pass == 0 && slice < ARGON2_SYNC_POINTS / 2;
#else
    data_independent = false;
#endif
    if (data_independent) {
        uint addr_index = offset % ARGON2_QWORDS_IN_BLOCK;
        if (addr_index == 0) {
            if (thread == 6) {
                ++*thread_input;
            }
            next_addresses(addr, tmp, *thread_input, thread, shuffle_buf);
        }

        uint thr = addr_index % THREADS_PER_LANE;
        uint idx = addr_index / THREADS_PER_LANE;

        ulong v = block_th_get(addr, idx);
        v = u64_shuffle(v, thr, thread, shuffle_buf);
        ref_index = u64_lo(v);
        ref_lane  = u64_hi(v);
    } else {
        ulong v = u64_shuffle(prev->a, 0, thread, shuffle_buf);
        ref_index = u64_lo(v);
        ref_lane  = u64_hi(v);
    }

    compute_ref_pos(lanes, segment_blocks, pass, lane, slice, offset,
                    &ref_lane, &ref_index);

    argon2_core(memory, mem_curr, prev, tmp, shuffle_buf, lanes, thread, pass,
                ref_index, ref_lane);
}

__kernel void argon2_kernel_segment(
        __local struct u64_shuffle_buf *shuffle_bufs,
        __global struct block_g *memory, uint passes, uint lanes,
        uint segment_blocks, uint pass, uint slice)
{
    uint job_id = get_global_id(1);
    uint lane   = get_global_id(0) / THREADS_PER_LANE;
    uint warp   = (get_local_id(1) * get_local_size(0) + get_local_id(0))
            / THREADS_PER_LANE;
    uint thread = get_local_id(0) % THREADS_PER_LANE;

    __local struct u64_shuffle_buf *shuffle_buf = &shuffle_bufs[warp];

    uint lane_blocks = ARGON2_SYNC_POINTS * segment_blocks;

    /* select job's memory region: */
    memory += (size_t)job_id * lanes * lane_blocks;

    struct block_th prev, addr, tmp;
    uint thread_input;

#if ARGON2_TYPE == ARGON2_I || ARGON2_TYPE == ARGON2_ID
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
        thread_input = ARGON2_TYPE;
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
#endif

    __global struct block_g *mem_segment =
            memory + slice * segment_blocks * lanes + lane;
    __global struct block_g *mem_prev, *mem_curr;
    uint start_offset = 0;
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

    for (uint offset = start_offset; offset < segment_blocks; ++offset) {
        argon2_step(memory, mem_curr, &prev, &tmp, &addr, shuffle_buf,
                    lanes, segment_blocks, thread, &thread_input,
                    lane, pass, slice, offset);

        mem_curr += lanes;
    }
}

__kernel void argon2_kernel_oneshot(
        __local struct u64_shuffle_buf *shuffle_bufs,
        __global struct block_g *memory, uint passes, uint lanes,
        uint segment_blocks)
{
    uint job_id = get_global_id(1);
    uint lane   = get_global_id(0) / THREADS_PER_LANE;
    uint warp   = get_local_id(1) * lanes + get_local_id(0) / THREADS_PER_LANE;
    uint thread = get_local_id(0) % THREADS_PER_LANE;

    __local struct u64_shuffle_buf *shuffle_buf = &shuffle_bufs[warp];

    uint lane_blocks = ARGON2_SYNC_POINTS * segment_blocks;

    /* select job's memory region: */
    memory += (size_t)job_id * lanes * lane_blocks;

    struct block_th prev, addr, tmp;
    uint thread_input;

#if ARGON2_TYPE == ARGON2_I || ARGON2_TYPE == ARGON2_ID
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
        thread_input = ARGON2_TYPE;
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
#endif

    __global struct block_g *mem_lane = memory + lane;
    __global struct block_g *mem_prev = mem_lane + 1 * lanes;
    __global struct block_g *mem_curr = mem_lane + 2 * lanes;

    load_block(&prev, mem_prev, thread);

    uint skip = 2;
    for (uint pass = 0; pass < passes; ++pass) {
        for (uint slice = 0; slice < ARGON2_SYNC_POINTS; ++slice) {
            for (uint offset = 0; offset < segment_blocks; ++offset) {
                if (skip > 0) {
                    --skip;
                    continue;
                }

                argon2_step(memory, mem_curr, &prev, &tmp, &addr, shuffle_buf,
                            lanes, segment_blocks, thread, &thread_input,
                            lane, pass, slice, offset);

                mem_curr += lanes;
            }

            barrier(CLK_GLOBAL_MEM_FENCE);

#if ARGON2_TYPE == ARGON2_I || ARGON2_TYPE == ARGON2_ID
            if (thread == 2) {
                ++thread_input;
            }
            if (thread == 6) {
                thread_input = 0;
            }
#endif
        }
#if ARGON2_TYPE == ARGON2_I
        if (thread == 0) {
            ++thread_input;
        }
        if (thread == 2) {
            thread_input = 0;
        }
#endif
        mem_curr = mem_lane;
    }
}

__kernel void argon2_kernel_preseed(
        __global struct block_g *memory, __global uint *seed, uint lanes, uint segment_blocks) {
    int job_id = get_global_id(1);
    int lane = get_global_id(0) % lanes;
    int idx = get_global_id(0) / lanes;

	/* select job's memory region: */
    memory += job_id * lanes * ARGON2_SYNC_POINTS * segment_blocks;
    seed += job_id * ARGON2_PREHASH_DIGEST_LENGTH;

    __global uint *initHash = (__global uint*)(memory + lane + (idx + 2) * lanes)->data;
    for(int i=0;i<ARGON2_PREHASH_DIGEST_LENGTH;i++) {
        initHash[i] = seed[i];
    }

    initHash[ARGON2_PREHASH_DIGEST_LENGTH] = idx;
    initHash[ARGON2_PREHASH_DIGEST_LENGTH + 1] = lane;
    blake2b_digestLong((__global uint*)(memory + lane + idx * lanes)->data, ARGON2_DWORDS_IN_BLOCK, initHash, ARGON2_PREHASH_SEED_LENGTH);
}

__kernel void argon2_kernel_finalize(
        __global struct block_g *memory, __global uint *out, uint outLen, uint lanes, uint segment_blocks) {
    int job_id = get_global_id(1);
    int thread = get_global_id(0);

    int lane_blocks = ARGON2_SYNC_POINTS * segment_blocks;
    /* select job's memory region: */
    memory += ((job_id + 1) * lanes * lane_blocks - lanes);
    out += job_id * outLen;
    __global struct block_g *dst = memory;

    for(int l=1;l<lanes;l++) {
        memory += 1;
        for (int i = 0; i < 4/*ARGON2_QWORDS_IN_BLOCK*/; i++) {
            dst->data[thread * 4 + i] ^= memory->data[thread * 4 + i];
        }
    }

    if(thread == 0) {
        blake2b_digestLong(out, outLen, (__global uint *) dst, ARGON2_DWORDS_IN_BLOCK);
    }
}
