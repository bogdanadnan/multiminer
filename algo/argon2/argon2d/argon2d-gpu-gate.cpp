#include "argon2d-gpu-gate.h"

#include "argon2-gpu/include/argon2-gpu-common/argon2params.h"
#include "argon2-gpu/include/argon2-opencl/processingunit.h"
#include "argon2-gpu/include/argon2-cuda/processingunit.h"
#include "argon2-gpu/include/argon2-cuda/cudaexception.h"

#include <iostream>
#include <strstream>
#include <unordered_map>
#include <mutex>
#include <vector>
#include <algorithm>

using namespace argon2;

int gpu_device_count = 0;
char *use_gpu = NULL;
std::vector<int> gpu_ids;
int gpu_batch_size = 1;
int total_threads = 1;

std::unordered_map<int, argon2_gpu_hasher_thread *> argon2_gpu_hashers;
std::mutex argon2_gpu_hashers_mutex;

std::vector<int> parse_gpu_id(const std::string &arg) {
	std::string::size_type pos, lastPos = 0, length = arg.length();
	std::vector<int> tokens;

	while(lastPos < length + 1)
	{
		pos = arg.find_first_of(",", lastPos);
		if(pos == std::string::npos)
		{
			pos = length;
		}

		if(pos != lastPos) {
			std::string token = std::string(arg.c_str() + lastPos,
								  pos - lastPos);
			int id = atoi(token.c_str()) - 1;
			if(id >= 0 && std::find(tokens.begin(), tokens.end(), id) == tokens.end())
				tokens.push_back(id);
		}

		lastPos = pos + 1;
	}

	std::sort(tokens.begin(), tokens.end());
	return tokens;
}

template<class GlobalContext, class ProgramContext, class ProcessingUnit>
bool init_gpu(int thr_id, CoinAlgo algo, Type type, Version version, Argon2Params *params) {
	argon2_gpu_hashers_mutex.lock();

	GlobalContext *context = new GlobalContext;
	auto &devices = context->getAllDevices();

	int selected_device = 0;
	if(gpu_ids.size() > 0) {
	    int selected_entry = thr_id / (total_threads / gpu_device_count);
	    if(selected_entry < 0 && selected_entry >= gpu_ids.size()) {
            argon2_gpu_hashers_mutex.unlock();
            return false;
        }
		selected_device = gpu_ids[selected_entry];
	}
	else {
		selected_device = thr_id / (total_threads / gpu_device_count);
	}

	for (std::size_t i = 0; i < devices.size(); i++) {
		if(i == selected_device) {
			auto &device = devices[i];

			argon2_gpu_hasher_thread *argon2_gpu_hasher_thread_data = (argon2_gpu_hasher_thread *)malloc(sizeof(argon2_gpu_hasher_thread));
			if(algo == Crds || algo == Dyn || algo == Arg || algo == Urx) {
				argon2_gpu_hasher_thread_data->vhash = (uint32_t *) malloc(gpu_batch_size * 8 * sizeof(uint32_t));
				argon2_gpu_hasher_thread_data->endiandata = (uint32_t *) malloc(20 * sizeof(uint32_t));
			}
			argon2_gpu_hasher_thread_data->processing_unit = NULL;

			ProgramContext *progCtx = new ProgramContext(context, {device}, type, version);
			argon2_gpu_hasher_thread_data->processing_unit = new ProcessingUnit(progCtx, params, &device, gpu_batch_size, algo, false);

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

		pu->setInput(thread_data->endiandata, 80);

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

bool init_thread_argon2d(int thr_id, argon2::CoinAlgo algo, argon2::Type type, argon2::Version version, Argon2Params *params) {
	if(use_gpu[0] == 'C') {
		try {
			if(!init_gpu<cuda::GlobalContext, cuda::ProgramContext, cuda::ProcessingUnit>(thr_id, algo, type, version, params)) {
				return false;
			}
		} catch (cuda::CudaException &err) {
			std::cerr << "CUDA ERROR: " << err.what() << std::endl;
			return false;
		}
	}
	else {
		try {
			if(!init_gpu<opencl::GlobalContext, opencl::ProgramContext, opencl::ProcessingUnit>(thr_id, algo, type, version, params)) {
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
	return init_thread_argon2d(thr_id, Arg, ARGON2_D, ARGON2_VERSION_13, new Argon2Params(32, nullptr, 0, nullptr, 0, nullptr, 0, 1, 4096, 1));
}

bool init_thread_argon2d_dyn_gpu(int thr_id) {
	return init_thread_argon2d(thr_id, Dyn, ARGON2_D, ARGON2_VERSION_10, new Argon2Params(32, nullptr, 0, nullptr, 0, nullptr, 0, 2, 500, 8));
}

bool init_thread_argon2d_crds_gpu(int thr_id) {
	return init_thread_argon2d(thr_id, Crds, ARGON2_D, ARGON2_VERSION_10, new Argon2Params(32, nullptr, 0, nullptr, 0, nullptr, 0, 1, 250, 4));
}

bool init_thread_argon2ad_urx_gpu_proxy(int thr_id, uint8_t *secret_ptr, size_t secret_sz, uint8_t *ad_ptr, size_t ad_sz) {
	return init_thread_argon2d(thr_id, Urx, ARGON2_D, ARGON2_VERSION_13, new Argon2Params(32, nullptr, 0, secret_ptr, secret_sz, ad_ptr, ad_sz, 1, 512, 2));
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

std::string join_ids(const std::vector<int> &ids) {
    std::ostrstream dest;
    for(int i=0; i < ids.size(); i++) {
        dest << "#" << (ids[i] + 1) << ((i < (ids.size() - 1)) ? ", " : "");
    }
    return dest.str();
}

template<class GlobalContext>
int show_gpu_info() {
	GlobalContext *context = new GlobalContext;
	auto &devices = context->getAllDevices();

	for (std::size_t i = 0; i < devices.size(); i++) {
		std::cout << "[Device #" << (i + 1) << "] " << devices[i].getName()
				  << std::endl;

	}

	bool gpu_id_set = !gpu_ids.empty();

	for(std::vector<int>::iterator it = gpu_ids.end(); it-- != gpu_ids.begin();) {
		if(*it < 0 || *it >= devices.size()) // invalid id, remove it
			gpu_ids.erase(it);
	}

	if(gpu_id_set && gpu_ids.size() == 0)
		std::cout << "Invalid GPU id passed in arguments, reverting to use all available devices." << std::endl;

	if(gpu_ids.size() == 0)
	    std::cout << "Start mining on all devices." << std::endl;
	else if(gpu_ids.size() == 1)
		std::cout << "Start mining on device #" << (gpu_ids[0] + 1) << "." << std::endl;
	else
        std::cout << "Start mining on devices " << join_ids(gpu_ids) << "." << std::endl;

	std::cout<<std::endl;

	gpu_device_count = gpu_ids.empty() ? devices.size() : gpu_ids.size();

	if(gpu_device_count > 0)
		total_threads *= gpu_device_count;

	return gpu_device_count;
}

int check_gpu_capability(char *_use_gpu, char *_gpu_id, int _gpu_batch_size, int threads) {
	if(_use_gpu == NULL) {
		return 0;
	}
	
	use_gpu = _use_gpu;

	if(_gpu_id != NULL) {
        gpu_ids = parse_gpu_id(_gpu_id);

        if (strlen(_gpu_id) > 0 && gpu_ids.size() == 0)
            std::cout << "Invalid GPU id passed in arguments, reverting to use all available devices." << std::endl;
    }

	gpu_batch_size = _gpu_batch_size;
	total_threads = threads;

	if(use_gpu[0] == 'C') {
		return show_gpu_info<cuda::GlobalContext>();
	}
	else {
		return show_gpu_info<opencl::GlobalContext>();
	}
}
