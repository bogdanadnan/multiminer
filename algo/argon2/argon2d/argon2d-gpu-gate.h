#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifndef WIN32
#define DLLEXPORT
#else
	#ifdef EXPORT_SYMBOLS
		#define DLLEXPORT __declspec(dllexport)
	#else
		#define DLLEXPORT __declspec(dllimport)
	#endif	
#endif	

typedef struct argon2_gpu_hasher_thread_ {
	void *processing_unit;
	uint32_t *vhash;
	uint32_t *endiandata;
} argon2_gpu_hasher_thread;

DLLEXPORT int check_gpu_capability(char *_use_gpu, int _gpu_id, int _gpu_batch_size, int threads);
DLLEXPORT argon2_gpu_hasher_thread *get_gpu_thread_data(int thr_id);
DLLEXPORT void gpu_argon2_raw_hash(argon2_gpu_hasher_thread *thread_data);
DLLEXPORT bool init_thread_argon2d4096_gpu(int thr_id);
DLLEXPORT bool init_thread_argon2d_dyn_gpu(int thr_id);
DLLEXPORT bool init_thread_argon2d_crds_gpu(int thr_id);

#ifdef __cplusplus
}
#endif
