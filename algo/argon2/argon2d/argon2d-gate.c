#include "argon2d-gate.h"
#include "argon2d/argon2.h"
#include <time.h>

static const size_t INPUT_BYTES = 80;  // Lenth of a block header in bytes. Input Length = Salt Length (salt = input)
static const size_t OUTPUT_BYTES = 32; // Length of output needed for a 256-bit hash
static const unsigned int DEFAULT_ARGON2_FLAG = 2; //Same as ARGON2_DEFAULT_FLAGS

bool init_thread_argon2d4096(int thr_id) {
	return init_thread_argon2d4096_gpu(thr_id);
}
bool init_thread_argon2d_dyn(int thr_id) {
	return init_thread_argon2d_dyn_gpu(thr_id);
}
bool init_thread_argon2d_crds(int thr_id) {
	return init_thread_argon2d_crds_gpu(thr_id);
}
// Credits

void argon2d_crds_hash( void *output, const void *input )
{
	argon2_context context;
	context.out = (uint8_t *)output;
	context.outlen = (uint32_t)OUTPUT_BYTES;
	context.pwd = (uint8_t *)input;
	context.pwdlen = (uint32_t)INPUT_BYTES;
	context.salt = (uint8_t *)input; //salt = input
	context.saltlen = (uint32_t)INPUT_BYTES;
	context.secret = NULL;
	context.secretlen = 0;
	context.ad = NULL;
	context.adlen = 0;
	context.allocate_cbk = NULL;
	context.free_cbk = NULL;
	context.flags = DEFAULT_ARGON2_FLAG; // = ARGON2_DEFAULT_FLAGS
	// main configurable Argon2 hash parameters
	context.m_cost = 250; // Memory in KiB (~256KB)
	context.lanes = 4;    // Degree of Parallelism
	context.threads = 1;  // Threads
	context.t_cost = 1;   // Iterations
        context.version = ARGON2_VERSION_10;

	argon2_ctx( &context, Argon2_d );
}

int scanhash_argon2d_crds( int thr_id, struct work *work, uint32_t max_nonce,
                      uint64_t *hashes_done )
{
        uint32_t _ALIGN(64) endiandata[20];
        uint32_t _ALIGN(64) hash[8];
        uint32_t *pdata = work->data;
        uint32_t *ptarget = work->target;

        const uint32_t first_nonce = pdata[19];
        const uint32_t Htarg = ptarget[7];

        uint32_t nonce = first_nonce;

        swab32_array( endiandata, pdata, 20 );

        do {
                be32enc(&endiandata[19], nonce);
                argon2d_crds_hash( hash, endiandata );
                if ( hash[7] <= Htarg && fulltest( hash, ptarget ) )
                {
                        pdata[19] = nonce;
                        *hashes_done = pdata[19] - first_nonce;
                        work_set_target_ratio(work, hash);
                        return 1;
                }
                nonce++;
        } while (nonce < max_nonce && !work_restart[thr_id].restart);

        pdata[19] = nonce;
        *hashes_done = pdata[19] - first_nonce + 1;
        return 0;
}

bool register_argon2d_crds_algo( algo_gate_t* gate )
{
	if(use_gpu == NULL) {
		gate->scanhash = (void *) &scanhash_argon2d_crds;
		gate->hash = (void *) &argon2d_crds_hash;
	}
	else {
		gate->miner_thread_pre_init = (void *) &init_thread_argon2d_crds;
		gate->scanhash = (void *) &scanhash_argon2d_crds_gpu;
	}
	gate->set_target = (void*)&scrypt_set_target;
	gate->optimizations = SSE2_OPT | AVX2_OPT | AVX512_OPT;
	return true;
}

// Dynamic

void argon2d_dyn_hash( void *output, const void *input )
{
    argon2_context context;
    context.out = (uint8_t *)output;
    context.outlen = (uint32_t)OUTPUT_BYTES;
    context.pwd = (uint8_t *)input;
    context.pwdlen = (uint32_t)INPUT_BYTES;
    context.salt = (uint8_t *)input; //salt = input
    context.saltlen = (uint32_t)INPUT_BYTES;
    context.secret = NULL;
    context.secretlen = 0;
    context.ad = NULL;
    context.adlen = 0;
    context.allocate_cbk = NULL;
    context.free_cbk = NULL;
    context.flags = DEFAULT_ARGON2_FLAG; // = ARGON2_DEFAULT_FLAGS
    // main configurable Argon2 hash parameters
    context.m_cost = 500;  // Memory in KiB (512KB)
    context.lanes = 8;     // Degree of Parallelism
    context.threads = 1;   // Threads
    context.t_cost = 2;    // Iterations
    context.version = ARGON2_VERSION_10;

    argon2_ctx( &context, Argon2_d );
}

int scanhash_argon2d_dyn( int thr_id, struct work *work, uint32_t max_nonce,
                      uint64_t *hashes_done )
{
        uint32_t _ALIGN(64) endiandata[20];
        uint32_t _ALIGN(64) hash[8];
        uint32_t *pdata = work->data;
        uint32_t *ptarget = work->target;

        const uint32_t first_nonce = pdata[19];
        const uint32_t Htarg = ptarget[7];

        uint32_t nonce = first_nonce;

        swab32_array( endiandata, pdata, 20 );

        do {
                be32enc(&endiandata[19], nonce);
                argon2d_dyn_hash( hash, endiandata );
                if ( hash[7] <= Htarg && fulltest( hash, ptarget ) )
                {
                        pdata[19] = nonce;
                        *hashes_done = pdata[19] - first_nonce;
                        work_set_target_ratio(work, hash);
                        return 1;
                }
                nonce++;
        } while (nonce < max_nonce && !work_restart[thr_id].restart);

        pdata[19] = nonce;
        *hashes_done = pdata[19] - first_nonce + 1;
        return 0;
}

bool register_argon2d_dyn_algo( algo_gate_t* gate )
{
		if(use_gpu == NULL) {
			gate->scanhash = (void *) &scanhash_argon2d_dyn;
			gate->hash = (void *) &argon2d_dyn_hash;
		}
		else {
			gate->miner_thread_pre_init = (void *) &init_thread_argon2d_dyn;
			gate->scanhash = (void *) &scanhash_argon2d_dyn_gpu;
		}
        gate->set_target = (void*)&scrypt_set_target;
        gate->optimizations = SSE2_OPT | AVX2_OPT | AVX512_OPT;
        return true;
}

// Unitus

int scanhash_argon2d4096( int thr_id, struct work *work, uint32_t max_nonce,
                           uint64_t *hashes_done)
{
   uint32_t _ALIGN(64) vhash[8];
   uint32_t _ALIGN(64) endiandata[20];
   uint32_t *pdata = work->data;
   uint32_t *ptarget = work->target;
   const uint32_t Htarg = ptarget[7];
   const uint32_t first_nonce = pdata[19];
   uint32_t n = first_nonce;
    
   uint32_t t_cost = 1; // 1 iteration
   uint32_t m_cost = 4096; // use 4MB
   uint32_t parallelism = 1; // 1 thread, 2 lanes

   for ( int i = 0; i < 19; i++ )
      be32enc( &endiandata[i], pdata[i] );

   do {
      be32enc( &endiandata[19], n );
      argon2d_hash_raw( t_cost, m_cost, parallelism, (char*) endiandata, 80,
                 (char*) endiandata, 80, (char*) vhash, 32, ARGON2_VERSION_13 );
      if ( vhash[7] < Htarg && fulltest( vhash, ptarget ) )
      {
         *hashes_done = n - first_nonce + 1;
         pdata[19] = n;
         return true;
      }
      n++;

   } while (n < max_nonce && !work_restart[thr_id].restart);

   *hashes_done = n - first_nonce + 1;
   pdata[19] = n;

   return 0;
}

int64_t get_max64_0x1ff() { return 0x1ff; }

bool register_argon2d4096_algo( algo_gate_t* gate )
{
    if(use_gpu == NULL)
        gate->scanhash = (void*)&scanhash_argon2d4096;
    else {
        gate->miner_thread_pre_init = (void *) &init_thread_argon2d4096;
        gate->scanhash = (void *) &scanhash_argon2d4096_gpu;
    }

    gate->set_target = (void*)&scrypt_set_target;
    gate->get_max64  = (void*)&get_max64_0x1ff;
    gate->optimizations = SSE2_OPT | AVX2_OPT | AVX512_OPT;
    return true;
}

// UraniumX

#define nullptr ((void*)0)

static const char* POW_SECRET = "f412a69fdc6d8ee6663f796b2e7ea53a52b9532a641b2f9cb2a7860108dc4c03";
static uint8_t* pArgon2Ad = nullptr;
static size_t nArgon2AdLen = 0;

static size_t Argon2FactorN(const int64_t nTime) {
	static const size_t offset = 9;
	static const int64_t nTimes[] = {
			0,          //                             512KB
			1618876800, // 04/20/2021 @ 12:00am (UTC)  1MB
			1713571200, // 04/20/2024 @ 12:00am (UTC)  2MB
			1808179200  // 04/20/2027 @ 12:00am (UTC)  4MB
	};
	size_t nFactor = 0;
	for (nFactor = 0; nFactor < 3; ++nFactor)
		if (nTime >= nTimes[nFactor] && nTime < nTimes[nFactor+1])
			return nFactor + offset;
	return nFactor + offset;
}

static uint32_t GetArgon2AdSize(const int64_t nTime) {
	int factor = Argon2FactorN(nTime);
	return 1024 * (1 << factor);
}

static void UpdateArgon2AdValues() {
	for (int i = 0; i < nArgon2AdLen; ++i)
		pArgon2Ad[i] = (uint8_t)(i < 256 ? i : i % 256);
}

static void EnsureArgon2MemoryAllocated(const int64_t nTime) {
	int nSize = GetArgon2AdSize(nTime);
	if (nSize > nArgon2AdLen) {
		if (nullptr != pArgon2Ad)
			free (pArgon2Ad);
		pArgon2Ad = (uint8_t*) malloc(nSize);
		nArgon2AdLen = nSize;
		UpdateArgon2AdValues();
	}
}

int Argon2Init() {
	int64_t nTime       = (int)time(NULL);
	EnsureArgon2MemoryAllocated(nTime);
	return (int)(nArgon2AdLen);
}

void Argon2Deinit() {
	if (nullptr != pArgon2Ad) {
		free(pArgon2Ad);
		pArgon2Ad = nullptr;
	}
}

void argon2ad_urx_hash( void *output, const void *input )
{
	argon2_context ctx;
	ctx.version         = ARGON2_VERSION_13;
	ctx.flags           = ARGON2_DEFAULT_FLAGS;
	ctx.out             = (uint8_t*) output;
	ctx.outlen          = OUTPUT_BYTES;
	ctx.pwd             = (uint8_t*)input;
	ctx.pwdlen          = INPUT_BYTES - 40;
	ctx.salt            = ((uint8_t*) input) + 40;
	ctx.saltlen         = 40;
	ctx.secret          = (uint8_t*) POW_SECRET;
	ctx.secretlen       = strlen (POW_SECRET);
	ctx.ad              = pArgon2Ad;
	ctx.adlen           = nArgon2AdLen;
	ctx.m_cost          = 512;
	ctx.t_cost          = 1;
	ctx.lanes           = 2;
	ctx.threads         = 1;
	ctx.allocate_cbk    = nullptr;
	ctx.free_cbk        = nullptr;
	argon2_ctx (&ctx, Argon2_d);
}

int scanhash_argon2ad_urx( int thr_id, struct work *work, uint32_t max_nonce,
						   uint64_t *hashes_done )
{
	uint32_t _ALIGN(64) endiandata[20];
	uint32_t _ALIGN(64) hash[8];
	uint32_t *pdata = work->data;
	uint32_t *ptarget = work->target;

	const uint32_t first_nonce = pdata[19];
	const uint32_t Htarg = ptarget[7];

	uint32_t nonce = first_nonce;

	swab32_array( endiandata, pdata, 20 );

	do {
		be32enc(&endiandata[19], nonce);
		argon2ad_urx_hash( hash, endiandata );
		if ( hash[7] <= Htarg && fulltest( hash, ptarget ) )
		{
			pdata[19] = nonce;
			*hashes_done = pdata[19] - first_nonce;
			work_set_target_ratio(work, hash);
			return 1;
		}
		nonce++;
	} while (nonce < max_nonce && !work_restart[thr_id].restart);

	pdata[19] = nonce;
	*hashes_done = pdata[19] - first_nonce + 1;
	return 0;
}

bool init_thread_argon2ad_urx_cpu(int thr_id) {
	return Argon2Init() > 0;
}

bool init_thread_argon2ad_urx_gpu(int thr_id) {
	bool success = (Argon2Init() > 0);
	if(success)
		success = init_thread_argon2ad_urx_gpu_proxy(thr_id, (uint8_t*) POW_SECRET, strlen (POW_SECRET), pArgon2Ad, nArgon2AdLen);

	return success;
}

bool register_argon2ad_urx_algo( algo_gate_t* gate )
{
	if(use_gpu == NULL) {
		gate->miner_thread_pre_init = (void *) &init_thread_argon2ad_urx_cpu;
		gate->scanhash = (void *) &scanhash_argon2ad_urx;
	}
	else {
		gate->miner_thread_pre_init = (void *) &init_thread_argon2ad_urx_gpu;
		gate->scanhash = (void *) &scanhash_argon2ad_urx_gpu;
	}
	gate->hash = (void*)&argon2ad_urx_hash;
	gate->set_target = (void*)&scrypt_set_target;
	gate->optimizations = SSE2_OPT | AVX2_OPT | AVX512_OPT;
	return true;
}

int scanhash_argon2d_dyn_gpu(int thr_id, struct work *work, uint32_t max_nonce,
							 uint64_t *hashes_done) {
	argon2_gpu_hasher_thread *thread_data = get_gpu_thread_data(thr_id);
	uint32_t *vhash;

	uint32_t *pdata = work->data;
	uint32_t *ptarget = work->target;
	const uint32_t Htarg = ptarget[7];
	const uint32_t first_nonce = pdata[19];
	uint32_t n = first_nonce;

	swab32_array(thread_data->endiandata, pdata, 20);

	do {
		thread_data->endiandata[19] = n;
		gpu_argon2_raw_hash(thread_data);

		for(int i=0;i<gpu_batch_size;i++) {
			vhash = thread_data->vhash + 8 * i;
			if (vhash[7] <= Htarg && fulltest(vhash, ptarget)) {
				*hashes_done = n - first_nonce;
				pdata[19] = n;
				work_set_target_ratio(work, vhash);
				return 1;
			}
			n++;
		}
	} while (n < max_nonce && !work_restart[thr_id].restart);

	*hashes_done = n - first_nonce + 1;
	pdata[19] = n;

	return 0;
}

int scanhash_argon2d4096_gpu(int thr_id, struct work *work, uint32_t max_nonce,
							 uint64_t *hashes_done) {
	argon2_gpu_hasher_thread *thread_data = get_gpu_thread_data(thr_id);
	uint32_t *vhash;

	uint32_t *pdata = work->data;
	uint32_t *ptarget = work->target;
	const uint32_t Htarg = ptarget[7];
	const uint32_t first_nonce = pdata[19];
	uint32_t n = first_nonce;

	swab32_array(thread_data->endiandata, pdata, 20);

	do {
		thread_data->endiandata[19] = n;
		gpu_argon2_raw_hash(thread_data);

		for(int i=0;i<gpu_batch_size;i++) {
			vhash = thread_data->vhash + 8 * i;
			if (vhash[7] < Htarg && fulltest(vhash, ptarget)) {
				*hashes_done = n - first_nonce + 1;
				pdata[19] = n;
				return 1;
			}
			n++;
		}
	} while (n < max_nonce && !work_restart[thr_id].restart);

	*hashes_done = n - first_nonce + 1;
	pdata[19] = n;

	return 0;
}

int scanhash_argon2d_crds_gpu(int thr_id, struct work *work, uint32_t max_nonce,
							 uint64_t *hashes_done) {
	argon2_gpu_hasher_thread *thread_data = get_gpu_thread_data(thr_id);
	uint32_t *vhash;

	uint32_t *pdata = work->data;
	uint32_t *ptarget = work->target;
	const uint32_t Htarg = ptarget[7];
	const uint32_t first_nonce = pdata[19];
	uint32_t n = first_nonce;

	swab32_array(thread_data->endiandata, pdata, 20);

	do {
		thread_data->endiandata[19] = n;
		gpu_argon2_raw_hash(thread_data);

		for(int i=0;i<gpu_batch_size;i++) {
			vhash = thread_data->vhash + 8 * i;
			if (vhash[7] <= Htarg && fulltest(vhash, ptarget)) {
				*hashes_done = n - first_nonce;
				pdata[19] = n;
				work_set_target_ratio(work, vhash);
				return 1;
			}
			n++;
		}
	} while (n < max_nonce && !work_restart[thr_id].restart);

	*hashes_done = n - first_nonce + 1;
	pdata[19] = n;

	return 0;
}

int scanhash_argon2ad_urx_gpu(int thr_id, struct work *work, uint32_t max_nonce,
							  uint64_t *hashes_done) {

	argon2_gpu_hasher_thread *thread_data = get_gpu_thread_data(thr_id);
	uint32_t *vhash;

	uint32_t *pdata = work->data;
	uint32_t *ptarget = work->target;
	const uint32_t Htarg = ptarget[7];
	const uint32_t first_nonce = pdata[19];
	uint32_t n = first_nonce;

	swab32_array(thread_data->endiandata, pdata, 20);

	do {
		thread_data->endiandata[19] = n;
		gpu_argon2_raw_hash(thread_data);

		for(int i=0;i<gpu_batch_size;i++) {
			vhash = thread_data->vhash + 8 * i;
			if (vhash[7] <= Htarg && fulltest(vhash, ptarget)) {
				*hashes_done = n - first_nonce;
				pdata[19] = n;
				work_set_target_ratio(work, vhash);
				return 1;
			}
			n++;
		}
	} while (n < max_nonce && !work_restart[thr_id].restart);

	*hashes_done = n - first_nonce + 1;
	pdata[19] = n;

	return 0;
}
