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

#define ROUND(m, t, r, shfl) \
do { \
	G(m, r, t, v0, v1, v2, v3); \
    shfl[t + 4] = v1; \
    shfl[t + 8] = v2; \
    shfl[t + 12] = v3; \
    mem_fence(CLK_LOCAL_MEM_FENCE); \
    v1 = shfl[((t + 1) % 4)+ 4]; \
    v2 = shfl[((t + 2) % 4)+ 8]; \
    v3 = shfl[((t + 3) % 4)+ 12]; \
	G(m, r, (t + 4), v0, v1, v2, v3); \
    shfl[((t + 1) % 4)+ 4] = v1; \
    shfl[((t + 2) % 4)+ 8] = v2; \
    shfl[((t + 3) % 4)+ 12] = v3; \
    mem_fence(CLK_LOCAL_MEM_FENCE); \
    v1 = shfl[t + 4]; \
    v2 = shfl[t + 8]; \
    v3 = shfl[t + 12]; \
} while ((void)0, 0)

ulong rotr64(ulong x, ulong n)
{
	return rotate(x, 64 - n);
}

__constant ulong blake2b_IV[8] = {
        0x6A09E667F3BCC908, 0xBB67AE8584CAA73B,
        0x3C6EF372FE94F82B, 0xA54FF53A5F1D36F1,
        0x510E527FADE682D1, 0x9B05688C2B3E6C1F,
        0x1F83D9ABFB41BD6B, 0x5BE0CD19137E2179
};

__constant uint blake2b_sigma[12][16] = {
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

void blake2b_compress(__local ulong *h, __local ulong *m, ulong f0, __local ulong *shfl, int thr_id)
{
    ulong v0, v1, v2, v3;

    mem_fence(CLK_LOCAL_MEM_FENCE);

    v0 = h[thr_id];
    v1 = h[thr_id + 4];
    v2 = blake2b_IV[thr_id];
    v3 = blake2b_IV[thr_id + 4];

    if(thr_id == 0) v3 ^= h[8];
    if(thr_id == 1) v3 ^= h[9];
    if(thr_id == 2) v3 ^= f0;

    ROUND(m, thr_id, 0, shfl);
    ROUND(m, thr_id, 1, shfl);
    ROUND(m, thr_id, 2, shfl);
    ROUND(m, thr_id, 3, shfl);
    ROUND(m, thr_id, 4, shfl);
    ROUND(m, thr_id, 5, shfl);
    ROUND(m, thr_id, 6, shfl);
    ROUND(m, thr_id, 7, shfl);
    ROUND(m, thr_id, 8, shfl);
    ROUND(m, thr_id, 9, shfl);
    ROUND(m, thr_id, 10, shfl);
    ROUND(m, thr_id, 11, shfl);

    h[thr_id] ^= v0 ^ v2;
    h[thr_id + 4] ^= v1 ^ v3;
}

void blake2b_incrementCounter(__local ulong *h, int inc)
{
    h[8] += (inc * 4);
    h[9] += (h[8] < (inc * 4));
}

void blake2b_final_global(__global uint *out, int out_len, __local ulong *h, __local uint *buf, int buf_len, __local ulong *shfl, int thr_id)
{
    int left = BLOCK_BYTES - buf_len;
    __local uint *cursor_out_local = buf + buf_len;

    for(int i=0; i < (left >> 2); i++, cursor_out_local += 4) {
        cursor_out_local[thr_id] = 0;
    }

    if(thr_id == 0) {
        for (int i = 0; i < (left % 4); i++) {
            cursor_out_local[i] = 0;
        }
        blake2b_incrementCounter(h, buf_len);
    }

    blake2b_compress(h, (__local ulong *)buf, 0xFFFFFFFFFFFFFFFF, shfl, thr_id);

    __local uint *cursor_in = (__local uint *)h;
    __global uint *cursor_out_global = out;

    for(int i=0; i < (out_len >> 2); i++, cursor_in += 4, cursor_out_global += 4) {
        cursor_out_global[thr_id] = cursor_in[thr_id];
    }

    if(thr_id == 0) {
        for (int i = 0; i < (out_len % 4); i++) {
            cursor_out_global[i] = cursor_in[i];
        }
    }
}

void blake2b_final_local(__local uint *out, int out_len, __local ulong *h, __local uint *buf, int buf_len, __local ulong *shfl, int thr_id)
{
    int left = BLOCK_BYTES - buf_len;
    __local uint *cursor_out = buf + buf_len;

    for(int i=0; i < (left >> 2); i++, cursor_out += 4) {
        cursor_out[thr_id] = 0;
    }

    if(thr_id == 0) {
        for (int i = 0; i < (left % 4); i++) {
            cursor_out[i] = 0;
        }
        blake2b_incrementCounter(h, buf_len);
    }

    blake2b_compress(h, (__local ulong *)buf, 0xFFFFFFFFFFFFFFFF, shfl, thr_id);

    __local uint *cursor_in = (__local uint *)h;
    cursor_out = out;

    for(int i=0; i < (out_len >> 2); i++, cursor_in += 4, cursor_out += 4) {
        cursor_out[thr_id] = cursor_in[thr_id];
    }

    if(thr_id == 0) {
        for (int i = 0; i < (out_len % 4); i++) {
            cursor_out[i] = cursor_in[i];
        }
    }
}

int blake2b_update_global(__global uint *in, int in_len, __local ulong *h, __local uint *buf, int buf_len, __local ulong *shfl, int thr_id)
{
    __global uint *cursor_in = in;
    __local uint *cursor_out = buf + buf_len;

    if (buf_len + in_len > BLOCK_BYTES) {
        int left = BLOCK_BYTES - buf_len;

        for(int i=0; i < (left >> 2); i++, cursor_in += 4, cursor_out += 4) {
            cursor_out[thr_id] = cursor_in[thr_id];
        }

        if(thr_id == 0) {
            for (int i = 0; i < (left % 4); i++) {
                cursor_out[i] = cursor_in[i];
            }
            blake2b_incrementCounter(h, BLOCK_BYTES);
        }

        blake2b_compress(h, (__local ulong *)buf, 0, shfl, thr_id);

        buf_len = 0;

        in_len -= left;
        in += left;

        while (in_len > BLOCK_BYTES) {
            if(thr_id == 0)
                blake2b_incrementCounter(h, BLOCK_BYTES);

            cursor_in = in;
            cursor_out = buf;

            for(int i=0; i < (BLOCK_BYTES / 4); i++, cursor_in += 4, cursor_out += 4) {
                cursor_out[thr_id] = cursor_in[thr_id];
            }

            blake2b_compress(h, (__local ulong *)buf, 0, shfl, thr_id);

            in_len -= BLOCK_BYTES;
            in += BLOCK_BYTES;
        }
    }

    cursor_in = in;
    cursor_out = buf + buf_len;

    for(int i=0; i < (in_len >> 2); i++, cursor_in += 4, cursor_out += 4) {
        cursor_out[thr_id] = cursor_in[thr_id];
    }

    if(thr_id == 0) {
        for (int i = 0; i < (in_len % 4); i++) {
            cursor_out[i] = cursor_in[i];
        }
    }

    return buf_len + in_len;
}

int blake2b_update_local(__local uint *in, int in_len, __local ulong *h, __local uint *buf, int buf_len, __local ulong *shfl, int thr_id)
{
    __local uint *cursor_in = in;
    __local uint *cursor_out = buf + buf_len;

    if (buf_len + in_len > BLOCK_BYTES) {
        int left = BLOCK_BYTES - buf_len;

        for(int i=0; i < (left >> 2); i++, cursor_in += 4, cursor_out += 4) {
            cursor_out[thr_id] = cursor_in[thr_id];
        }

        if(thr_id == 0) {
            for (int i = 0; i < (left % 4); i++) {
                cursor_out[i] = cursor_in[i];
            }
            blake2b_incrementCounter(h, BLOCK_BYTES);
        }

        blake2b_compress(h, (__local ulong *)buf, 0, shfl, thr_id);

        buf_len = 0;

        in_len -= left;
        in += left;

        while (in_len > BLOCK_BYTES) {
            if(thr_id == 0)
                blake2b_incrementCounter(h, BLOCK_BYTES);

            cursor_in = in;
            cursor_out = buf;

            for(int i=0; i < (BLOCK_BYTES / 4); i++, cursor_in += 4, cursor_out += 4) {
                cursor_out[thr_id] = cursor_in[thr_id];
            }

            blake2b_compress(h, (__local ulong *)buf, 0, shfl, thr_id);

            in_len -= BLOCK_BYTES;
            in += BLOCK_BYTES;
        }
    }

    cursor_in = in;
    cursor_out = buf + buf_len;

    for(int i=0; i < (in_len >> 2); i++, cursor_in += 4, cursor_out += 4) {
        cursor_out[thr_id] = cursor_in[thr_id];
    }

    if(thr_id == 0) {
        for (int i = 0; i < (in_len % 4); i++) {
            cursor_out[i] = cursor_in[i];
        }
    }

    return buf_len + in_len;
}

int blake2b_init(__local ulong *h, int out_len, int thr_id)
{
    h[thr_id * 2] = blake2b_IV[thr_id * 2];
    h[thr_id * 2 + 1] = blake2b_IV[thr_id * 2 + 1];

    if(thr_id == 0) {
        h[8] = h[9] = 0;
        h[0] = 0x6A09E667F3BCC908 ^ ((out_len * 4) | (1 << 16) | (1 << 24));
    }

    return 0;
}

void blake2b_digestLong_global(__global uint *out, int out_len,
                       __global uint *in, int in_len,
                       int thr_id, __local ulong* shared)
{
    __local ulong *h = shared;
	__local ulong *shfl = &h[10];
    __local uint *buf = (__local uint *)&shfl[16];
    __local uint *out_buffer = &buf[32];
    int buf_len;

    if(thr_id == 0) buf[0] = (out_len * 4);
    buf_len = 1;

    if (out_len <= OUT_BYTES) {
        blake2b_init(h, out_len, thr_id);
        buf_len = blake2b_update_global(in, in_len, h, buf, buf_len, shfl, thr_id);
        blake2b_final_global(out, out_len, h, buf, buf_len, shfl, thr_id);
    } else {
        __local uint *cursor_in = out_buffer;
        __global uint *cursor_out = out;

        blake2b_init(h, OUT_BYTES, thr_id);
        buf_len = blake2b_update_global(in, in_len, h, buf, buf_len, shfl, thr_id);
        blake2b_final_local(out_buffer, OUT_BYTES, h, buf, buf_len, shfl, thr_id);

        for(int i=0; i < (OUT_BYTES / 8); i++, cursor_in += 4, cursor_out += 4) {
            cursor_out[thr_id] = cursor_in[thr_id];
        }

        out += OUT_BYTES / 2;

        int to_produce = out_len - OUT_BYTES / 2;
        while (to_produce > OUT_BYTES) {
            buf_len = blake2b_init(h, OUT_BYTES, thr_id);
            buf_len = blake2b_update_local(out_buffer, OUT_BYTES, h, buf, buf_len, shfl, thr_id);
            blake2b_final_local(out_buffer, OUT_BYTES, h, buf, buf_len, shfl, thr_id);

            cursor_out = out;
            cursor_in = out_buffer;
            for(int i=0; i < (OUT_BYTES / 8); i++, cursor_in += 4, cursor_out += 4) {
                cursor_out[thr_id] = cursor_in[thr_id];
            }

            out += OUT_BYTES / 2;
            to_produce -= OUT_BYTES / 2;
        }

        buf_len = blake2b_init(h, to_produce, thr_id);
        buf_len = blake2b_update_local(out_buffer, OUT_BYTES, h, buf, buf_len, shfl, thr_id);
        blake2b_final_global(out, to_produce, h, buf, buf_len, shfl, thr_id);
    }
}

void blake2b_digestLong_local(__global uint *out, int out_len,
                        __local uint *in, int in_len,
                        int thr_id, __local ulong* shared)
{
    __local ulong *h = shared;
    __local ulong *shfl = &h[10];
    __local uint *buf = (__local uint *)&shfl[16];
    __local uint *out_buffer = &buf[32];
    int buf_len;

    if(thr_id == 0) buf[0] = (out_len * 4);
    buf_len = 1;

    if (out_len <= OUT_BYTES) {
        blake2b_init(h, out_len, thr_id);
        buf_len = blake2b_update_local(in, in_len, h, buf, buf_len, shfl, thr_id);
        blake2b_final_global(out, out_len, h, buf, buf_len, shfl, thr_id);
    } else {
        __local uint *cursor_in = out_buffer;
        __global uint *cursor_out = out;

        blake2b_init(h, OUT_BYTES, thr_id);
        buf_len = blake2b_update_local(in, in_len, h, buf, buf_len, shfl, thr_id);
        blake2b_final_local(out_buffer, OUT_BYTES, h, buf, buf_len, shfl, thr_id);

        for(int i=0; i < (OUT_BYTES / 8); i++, cursor_in += 4, cursor_out += 4) {
            cursor_out[thr_id] = cursor_in[thr_id];
        }

        out += OUT_BYTES / 2;

        int to_produce = out_len - OUT_BYTES / 2;
        while (to_produce > OUT_BYTES) {
            buf_len = blake2b_init(h, OUT_BYTES, thr_id);
            buf_len = blake2b_update_local(out_buffer, OUT_BYTES, h, buf, buf_len, shfl, thr_id);
            blake2b_final_local(out_buffer, OUT_BYTES, h, buf, buf_len, shfl, thr_id);

            cursor_out = out;
            cursor_in = out_buffer;
            for(int i=0; i < (OUT_BYTES / 8); i++, cursor_in += 4, cursor_out += 4) {
                cursor_out[thr_id] = cursor_in[thr_id];
            }

            out += OUT_BYTES / 2;
            to_produce -= OUT_BYTES / 2;
        }

        buf_len = blake2b_init(h, to_produce, thr_id);
        buf_len = blake2b_update_local(out_buffer, OUT_BYTES, h, buf, buf_len, shfl, thr_id);
        blake2b_final_global(out, to_produce, h, buf, buf_len, shfl, thr_id);
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

void argon2_genseed_generic(__local uint *initHash, __global uint *seed, int job_id, int thr_id) {
    __global uint *seed_local = seed + job_id * ARGON2_PREHASH_DIGEST_LENGTH;

    for (int i = 0; i < ARGON2_PREHASH_DIGEST_LENGTH / 4; i++) {
        initHash[i * 4 + thr_id] = seed_local[i * 4 + thr_id];
    }
}

void argon2_genseed_crds_dyn_arg(__local uint *initHash, __global uint *seed,
                                 int lanes, int m_cost, int t_cost, int version, int job_id, int thr_id) {
    __local ulong *h = (__local ulong *)&initHash[20];
    __local ulong *shfl = &h[10];
    __local uint *buf = (__local uint *)&shfl[16];
    __local uint *value = &buf[32];

    for (int i = 0; i < 5; i++) {
        initHash[i * 4 + thr_id] = seed[i * 4 + thr_id];
    }

    if (thr_id == 3) {
        uint x = seed[19] + job_id;
        __local uchar *p = (__local uchar *)&initHash[19];
        p[3] = x & 0xff;
        p[2] = (x >> 8) & 0xff;
        p[1] = (x >> 16) & 0xff;
        p[0] = (x >> 24) & 0xff;
    }

    int buf_len = blake2b_init(h, ARGON2_PREHASH_DIGEST_LENGTH, thr_id);
    *value = lanes; //lanes
    buf_len = blake2b_update_local(value, 1, h, buf, buf_len, shfl, thr_id);
    *value = 32; //outlen
    buf_len = blake2b_update_local(value, 1, h, buf, buf_len, shfl, thr_id);
    *value = m_cost; //m_cost
    buf_len = blake2b_update_local(value, 1, h, buf, buf_len, shfl, thr_id);
    *value = t_cost; //t_cost
    buf_len = blake2b_update_local(value, 1, h, buf, buf_len, shfl, thr_id);
    *value = version; //version
    buf_len = blake2b_update_local(value, 1, h, buf, buf_len, shfl, thr_id);
    *value = ARGON2_D; //type
    buf_len = blake2b_update_local(value, 1, h, buf, buf_len, shfl, thr_id);
    *value = 80; //pw_len
    buf_len = blake2b_update_local(value, 1, h, buf, buf_len, shfl, thr_id);
    buf_len = blake2b_update_local(initHash, 20, h, buf, buf_len, shfl, thr_id);
    *value = 80; //salt_len
    buf_len = blake2b_update_local(value, 1, h, buf, buf_len, shfl, thr_id);
    buf_len = blake2b_update_local(initHash, 20, h, buf, buf_len, shfl, thr_id);
    *value = 0; //secret_len
    buf_len = blake2b_update_local(value, 1, h, buf, buf_len, shfl, thr_id);
    buf_len = blake2b_update_local(0, 0, h, buf, buf_len, shfl, thr_id);
    *value = 0; //ad_len
    buf_len = blake2b_update_local(value, 1, h, buf, buf_len, shfl, thr_id);
    buf_len = blake2b_update_local(0, 0, h, buf, buf_len, shfl, thr_id);

    blake2b_final_local(initHash, ARGON2_PREHASH_DIGEST_LENGTH, h, buf, buf_len, shfl, thr_id);
}

void argon2_genseed_urx(__local uint *initHash, __global uint *seed, __global uint *secret, uint secretLen, __global uint *ad, uint adLen,
								 int lanes, int m_cost, int t_cost, int version, int job_id, int thr_id) {
	__local ulong *h = (__local ulong *)&initHash[20];
	__local ulong *shfl = &h[10];
	__local uint *buf = (__local uint *)&shfl[16];
	__local uint *value = &buf[32];

	for (int i = 0; i < 5; i++) {
		initHash[i * 4 + thr_id] = seed[i * 4 + thr_id];
	}

	if (thr_id == 3) {
		uint x = seed[19] + job_id;
		__local uchar *p = (__local uchar *)&initHash[19];
		p[3] = x & 0xff;
		p[2] = (x >> 8) & 0xff;
		p[1] = (x >> 16) & 0xff;
		p[0] = (x >> 24) & 0xff;
	}

	int buf_len = blake2b_init(h, ARGON2_PREHASH_DIGEST_LENGTH, thr_id);
	*value = lanes; //lanes
	buf_len = blake2b_update_local(value, 1, h, buf, buf_len, shfl, thr_id);
	*value = 32; //outlen
	buf_len = blake2b_update_local(value, 1, h, buf, buf_len, shfl, thr_id);
	*value = m_cost; //m_cost
	buf_len = blake2b_update_local(value, 1, h, buf, buf_len, shfl, thr_id);
	*value = t_cost; //t_cost
	buf_len = blake2b_update_local(value, 1, h, buf, buf_len, shfl, thr_id);
	*value = version; //version
	buf_len = blake2b_update_local(value, 1, h, buf, buf_len, shfl, thr_id);
	*value = ARGON2_D; //type
	buf_len = blake2b_update_local(value, 1, h, buf, buf_len, shfl, thr_id);
	*value = 40; //pw_len
	buf_len = blake2b_update_local(value, 1, h, buf, buf_len, shfl, thr_id);
	buf_len = blake2b_update_local(initHash, 10, h, buf, buf_len, shfl, thr_id);
	*value = 40; //salt_len
	buf_len = blake2b_update_local(value, 1, h, buf, buf_len, shfl, thr_id);
	buf_len = blake2b_update_local(&initHash[10], 10, h, buf, buf_len, shfl, thr_id);
	*value = secretLen; //secret_len
	buf_len = blake2b_update_local(value, 1, h, buf, buf_len, shfl, thr_id);
	buf_len = blake2b_update_global(secret, secretLen / 4, h, buf, buf_len, shfl, thr_id);
	*value = adLen; //ad_len
	buf_len = blake2b_update_local(value, 1, h, buf, buf_len, shfl, thr_id);
	buf_len = blake2b_update_global(ad, adLen / 4, h, buf, buf_len, shfl, thr_id);

	blake2b_final_local(initHash, ARGON2_PREHASH_DIGEST_LENGTH, h, buf, buf_len, shfl, thr_id);
}

__kernel void argon2_kernel_preseed(uint algo,
        __global struct block_g *memory, __global uint *seed, uint lanes, uint segment_blocks,
		__global uint *secret, uint secretLen, __global uint *ad, uint adLen,
		__local ulong *blake_shared) {
    int job_id = get_global_id(1);
    int lane_thr = get_global_id(0) / 4;
    int thr_id = get_global_id(0) % 4;
    int lane = lane_thr % lanes;
    int idx = lane_thr / lanes;

    __local uint *initHash = (__local uint *)&blake_shared[lane_thr * 60];

    if(algo == 1) // Crds
        argon2_genseed_crds_dyn_arg(initHash, seed, lanes, 250, 1, ARGON2_VERSION_10, job_id, thr_id);
    else if(algo == 2) // Dyn
        argon2_genseed_crds_dyn_arg(initHash, seed, lanes, 500, 2, ARGON2_VERSION_10, job_id, thr_id);
    else if(algo == 3) // Arg
        argon2_genseed_crds_dyn_arg(initHash, seed, lanes, 4096, 1, ARGON2_VERSION_13, job_id, thr_id);
    else if(algo == 4) //Urx
		argon2_genseed_urx(initHash, seed, secret, secretLen, ad, adLen, lanes, 512, 1, ARGON2_VERSION_13, job_id, thr_id);
    else
        argon2_genseed_generic(initHash, seed, job_id, thr_id);

    if (thr_id == 0) {
        initHash[ARGON2_PREHASH_DIGEST_LENGTH] = idx;
        initHash[ARGON2_PREHASH_DIGEST_LENGTH + 1] = lane;
    }

	/* select job's memory region: */
    memory += job_id * lanes * ARGON2_SYNC_POINTS * segment_blocks;

    blake2b_digestLong_local((__global uint*)(memory + lane + idx * lanes)->data, ARGON2_DWORDS_IN_BLOCK, initHash,
            ARGON2_PREHASH_SEED_LENGTH, thr_id, (__local ulong *)&initHash[20]);
}

__kernel void argon2_kernel_finalize(
        __global struct block_g *memory, __global uint *out, uint outLen, uint lanes, uint segment_blocks, __local ulong *blake_shared) {
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

    if(thread / 4 == 0) {
        blake2b_digestLong_global(out, outLen, (__global uint *) dst, ARGON2_DWORDS_IN_BLOCK, thread, blake_shared);
    }
}
