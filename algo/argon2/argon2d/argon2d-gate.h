#ifndef ARGON2D_GATE_H__
#define ARGON2D_GATE_H__

#include "algo-gate-api.h"
#include "argon2d-gpu-gate.h"

// Alterdot: version = 0x10, m_cost = 16000.
bool register_argon2d16000_algo( algo_gate_t* gate );

void argon2d16000_hash( void *state, const void *input );

int scanhash_argon2d16000( int thr_id, struct work *work, uint32_t max_nonce,
                    uint64_t *hashes_done );

// Credits / Zumy: version = 0x10, m_cost = 250.
bool register_argon2d_crds_algo( algo_gate_t* gate );

void argon2d_crds_hash( void *state, const void *input );

int scanhash_argon2d_crds( int thr_id, struct work *work, uint32_t max_nonce,
                    uint64_t *hashes_done );

// Dynamic: version = 0x10, m_cost = 500.
bool register_argon2d_dyn_algo( algo_gate_t* gate );

void argon2d_dyn_hash( void *state, const void *input );

int scanhash_argon2d_dyn( int thr_id, struct work *work, uint32_t max_nonce,
                    uint64_t *hashes_done );


// Unitus: version = 0x13, m_cost = 4096.
bool register_argon2d4096_algo( algo_gate_t* gate );

int scanhash_argon2d4096( int thr_id, struct work *work, uint32_t max_nonce,
						  uint64_t *hashes_done );

// functions to hash on GPU

// Unitus: version = 0x13, m_cost = 4096.
bool init_thread_argon2d4096(int thr_id);
int scanhash_argon2d4096_gpu(int thr_id, struct work *work, uint32_t max_nonce,
							  uint64_t *hashes_done);

// Dynamic: version = 0x10, m_cost = 500.
bool init_thread_argon2d_dyn(int thr_id);

// Credits / Zumy: version = 0x10, m_cost = 250.
bool init_thread_argon2d_crds(int thr_id);

// Alterdot: version = 0x10, m_cost = 16000.
bool init_thread_argon2d16000(int thr_id);

int scanhash_argon2d_multi_gpu( int thr_id, struct work *work, uint32_t max_nonce,
							  uint64_t *hashes_done );

#endif

