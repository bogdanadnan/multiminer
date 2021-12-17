#include "argon2d-gpu-gate.h"

#include "argon2-gpu/include/argon2-gpu-common/argon2params.h"
#include "argon2-gpu/include/argon2-opencl/processingunit.h"
#include "argon2-gpu/include/argon2-cuda/processingunit.h"
#include "argon2-gpu/include/argon2-cuda/cudaexception.h"

#include <iostream>
#include <unordered_map>
#include <mutex>

using namespace argon2;

int gpu_device_count = 0;
char *use_gpu = NULL;
int gpu_id = -1;
int gpu_batch_size = 1;
int total_threads = 1;

std::unordered_map<int, argon2_gpu_hasher_thread *> argon2_gpu_hashers;
std::mutex argon2_gpu_hashers_mutex;

template<class GlobalContext, class ProgramContext, class ProcessingUnit>
bool init_gpu(int thr_id, Type type, Version version, Argon2Params *params) {
	argon2_gpu_hashers_mutex.lock();

	GlobalContext *context = new GlobalContext;
	auto &devices = context->getAllDevices();

	int selected_device = 0;
	if(gpu_id >= 0) {
		selected_device = gpu_id;
	}
	else {
		selected_device = thr_id / (total_threads / gpu_device_count);
	}

	for (std::size_t i = 0; i < devices.size(); i++) {
		if(i == selected_device) {
			auto &device = devices[i];

			argon2_gpu_hasher_thread *argon2_gpu_hasher_thread_data = (argon2_gpu_hasher_thread *)malloc(sizeof(argon2_gpu_hasher_thread));
			argon2_gpu_hasher_thread_data->vhash = (uint32_t*)malloc(gpu_batch_size * 8 * sizeof(uint32_t));
			argon2_gpu_hasher_thread_data->endiandata = (uint32_t*)malloc(gpu_batch_size * 20 * sizeof(uint32_t));
			argon2_gpu_hasher_thread_data->processing_unit = NULL;

			ProgramContext *progCtx = new ProgramContext(context, {device}, type, version);
			argon2_gpu_hasher_thread_data->processing_unit = new ProcessingUnit(progCtx, params, &device, gpu_batch_size, false);

			std::cout << "[Thread " << thr_id << "] Device #" << (i + 1) << ": " << device.getName()
					  << std::endl << std::flush;

			argon2_gpu_hashers[thr_id] = argon2_gpu_hasher_thread_data;
		}
	}

	argon2_gpu_hashers_mutex.unlock();
	return true;
}

template<class ProcessingUnit>
void gpu_argon2_raw_hash_gate(argon2_gpu_hasher_thread *thread_data) {
	if(thread_data != NULL && thread_data->processing_unit != NULL) {
		ProcessingUnit *pu = (ProcessingUnit *) thread_data->processing_unit;

		for (std::size_t i = 0; i < gpu_batch_size; i++) {
			pu->setPasswordSameSalt(i, thread_data->endiandata + 20 * i, 80);
		}

		pu->beginProcessing();
		pu->endProcessing();

		for (std::size_t i = 0; i < gpu_batch_size; i++) {
			void *hash = pu->getHash(i);
			if(hash != NULL)
				memcpy(thread_data->vhash + 8 * i, hash, pu->getParams()->getOutputLength());
		}
	}
	else {
		std::cerr<<"Thread not initialized."<<std::endl;
		exit(1);
	}
}

void gpu_argon2_raw_hash(argon2_gpu_hasher_thread *thread_data) {
	if(use_gpu[0] == 'C')
		gpu_argon2_raw_hash_gate<cuda::ProcessingUnit>(thread_data);
	else
		gpu_argon2_raw_hash_gate<opencl::ProcessingUnit>(thread_data);
}

bool init_thread_argon2d(int thr_id, argon2::Type type, argon2::Version version, Argon2Params *params) {
	if(use_gpu[0] == 'C') {
		try {
			if(!init_gpu<cuda::GlobalContext, cuda::ProgramContext, cuda::ProcessingUnit>(thr_id, type, version, params)) {
				return false;
			}
		} catch (cuda::CudaException &err) {
			std::cerr << "CUDA ERROR: " << err.what() << std::endl;
			return false;
		}
	}
	else {
		try {
			if(!init_gpu<opencl::GlobalContext, opencl::ProgramContext, opencl::ProcessingUnit>(thr_id, type, version, params)) {
				return false;
			}
		} catch (cl::Error &err) {
			std::cerr << ": OpenCL ERROR: " << err.err() << ": "
					  << err.what() << std::endl;
			return false;
		}
	}
	return true;
}

bool init_thread_argon2d4096_gpu(int thr_id) {
	return init_thread_argon2d(thr_id, ARGON2_D, ARGON2_VERSION_13, new Argon2Params(32, nullptr, 0, nullptr, 0, nullptr, 0, 1, 4096, 1));
}

bool init_thread_argon2d_dyn_gpu(int thr_id) {
	return init_thread_argon2d(thr_id, ARGON2_D, ARGON2_VERSION_10, new Argon2Params(32, nullptr, 0, nullptr, 0, nullptr, 0, 2, 500, 8));
}

bool init_thread_argon2d_crds_gpu(int thr_id) {
	return init_thread_argon2d(thr_id, ARGON2_D, ARGON2_VERSION_10, new Argon2Params(32, nullptr, 0, nullptr, 0, nullptr, 0, 1, 250, 4));
}

bool init_thread_argon2d16000_gpu(int thr_id) {
	return init_thread_argon2d(thr_id, ARGON2_D, ARGON2_VERSION_10, new Argon2Params(32, nullptr, 0, nullptr, 0, nullptr, 0, 1, 16000, 1));
}

argon2_gpu_hasher_thread *get_gpu_thread_data(int thr_id) {
	argon2_gpu_hasher_thread *thread_data = NULL;
	argon2_gpu_hashers_mutex.lock();
	if(argon2_gpu_hashers.find(thr_id) != argon2_gpu_hashers.end()) {
		thread_data = argon2_gpu_hashers[thr_id];
	}
	argon2_gpu_hashers_mutex.unlock();
	return thread_data;
}

template<class GlobalContext>
int show_gpu_info() {
	GlobalContext *context = new GlobalContext;
	auto &devices = context->getAllDevices();

	for (std::size_t i = 0; i < devices.size(); i++) {
		std::cout << "[Device #" << (i + 1) << "] " << devices[i].getName()
				  << std::endl;

	}
	if(gpu_id >= 0)
		std::cout << "Start mining on device #" << (gpu_id + 1) << "." << std::endl;
	else
		std::cout << "Start mining on all devices." << std::endl;

	std::cout<<std::endl;

	gpu_device_count = devices.size();
	if(gpu_device_count > 0)
		total_threads *= gpu_device_count;

	return gpu_device_count;
}

int check_gpu_capability(char *_use_gpu, int _gpu_id, int _gpu_batch_size, int threads) {
	if(_use_gpu == NULL) {
		return 0;
	}
	
	use_gpu = _use_gpu;
	gpu_id = _gpu_id;
	gpu_batch_size = _gpu_batch_size;
	total_threads = threads;

	if(use_gpu[0] == 'C') {
		return show_gpu_info<cuda::GlobalContext>();
	}
	else {
		return show_gpu_info<opencl::GlobalContext>();
	}
}
