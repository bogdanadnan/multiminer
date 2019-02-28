#ifdef _WIN64
#include <windows.h>
#endif

#include "argon2-gpu-common/argon2-common.h"
#include "kernelrunner.h"

#include <stdexcept>

#ifndef NDEBUG
#include <iostream>
#endif

#define THREADS_PER_LANE 32

namespace argon2 {
namespace opencl {

enum {
    ARGON2_REFS_PER_BLOCK = ARGON2_BLOCK_SIZE / (2 * sizeof(cl_uint)),
};

KernelRunner::KernelRunner(const ProgramContext *programContext,
                           const Argon2Params *params, const Device *device,
                           std::size_t batchSize, size_t outLen,
                           bool bySegment, bool precompute,
						   std::uint8_t *secret_, std::size_t secretLen_,
						   std::uint8_t *ad_, std::size_t adLen_)
    : programContext(programContext), params(params), batchSize(batchSize), outLen(outLen),
      bySegment(bySegment), precompute(precompute),
      memorySize(params->getMemorySize() * batchSize)
{
    auto context = programContext->getContext();
    std::uint32_t passes = params->getTimeCost();
    std::uint32_t lanes = params->getLanes();
    std::uint32_t segmentBlocks = params->getSegmentBlocks();

    queue = cl::CommandQueue(context, device->getCLDevice(),
                             CL_QUEUE_PROFILING_ENABLE);

#ifndef NDEBUG
        std::cerr << "[INFO] Allocating " << memorySize << " bytes for memory..."
                  << std::endl;
#endif

    memoryBuffer = cl::Buffer(context, CL_MEM_READ_WRITE, memorySize);
	seedBuffer = cl::Buffer(context,CL_MEM_READ_ONLY, batchSize * ARGON2_PREHASH_DIGEST_LENGTH);
	outBuffer = cl::Buffer(context,CL_MEM_WRITE_ONLY, batchSize * outLen);
	seedHost = NULL;
	outHost = NULL;

    if(secret_ != NULL) {
        secretLen = secretLen_;
        secretBuffer = cl::Buffer(context, CL_MEM_READ_ONLY, secretLen);
        queue.enqueueWriteBuffer(secretBuffer, true, 0, secretLen, secret_);
    }
    else {
        secretLen = 0;
        secretBuffer = cl::Buffer(context, CL_MEM_READ_ONLY, 1); // dummy buffer
    }

    if(ad_ != NULL) {
        adLen = adLen_;
        adBuffer = cl::Buffer(context, CL_MEM_READ_ONLY, adLen);
        queue.enqueueWriteBuffer(adBuffer, true, 0, adLen, ad_);
    }
    else {
        adLen = 0;
        adBuffer = cl::Buffer(context, CL_MEM_READ_ONLY, 1); // dummy buffer
    }

    Type type = programContext->getArgon2Type();
    if ((type == ARGON2_I || type == ARGON2_ID) && precompute) {
        std::uint32_t segments =
                type == ARGON2_ID
                ? lanes * (ARGON2_SYNC_POINTS / 2)
                : passes * lanes * ARGON2_SYNC_POINTS;

        std::size_t refsSize = segments * segmentBlocks * sizeof(cl_uint) * 2;

#ifndef NDEBUG
        std::cerr << "[INFO] Allocating " << refsSize << " bytes for refs..."
                  << std::endl;
#endif

        refsBuffer = cl::Buffer(context, CL_MEM_READ_WRITE, refsSize);

        precomputeRefs();
    }

    static const char *KERNEL_NAMES[2][4] = {
        {
            "argon2_kernel_oneshot",
            "argon2_kernel_segment",
        },
        {
            "argon2_kernel_oneshot_precompute",
            "argon2_kernel_segment_precompute",
        }
    };

    kpreseed = cl::Kernel(programContext->getProgram(),
                         "argon2_kernel_preseed");

    kworker = cl::Kernel(programContext->getProgram(),
                         KERNEL_NAMES[precompute][bySegment]);

    kfinalize = cl::Kernel(programContext->getProgram(),
                         "argon2_kernel_finalize");

    kpreseed.setArg<cl::Buffer>(1, memoryBuffer);
    kpreseed.setArg<cl::Buffer>(2, seedBuffer);
    kpreseed.setArg<cl_uint>(3, lanes);
    kpreseed.setArg<cl_uint>(4, segmentBlocks);
    kpreseed.setArg<cl::Buffer>(5, secretBuffer);
    kpreseed.setArg<cl_uint>(6, secretLen);
    kpreseed.setArg<cl::Buffer>(7, adBuffer);
    kpreseed.setArg<cl_uint>(8, adLen);

    kworker.setArg<cl::Buffer>(1, memoryBuffer);
    if (precompute) {
        kworker.setArg<cl::Buffer>(2, refsBuffer);
        kworker.setArg<cl_uint>(3, passes);
        kworker.setArg<cl_uint>(4, lanes);
        kworker.setArg<cl_uint>(5, segmentBlocks);
    } else {
        kworker.setArg<cl_uint>(2, passes);
        kworker.setArg<cl_uint>(3, lanes);
        kworker.setArg<cl_uint>(4, segmentBlocks);
    }

    kfinalize.setArg<cl::Buffer>(0, memoryBuffer);
    kfinalize.setArg<cl::Buffer>(1, outBuffer);
    kfinalize.setArg<cl_uint>(2, outLen / 4);
    kfinalize.setArg<cl_uint>(3, lanes);
    kfinalize.setArg<cl_uint>(4, segmentBlocks);
}

void KernelRunner::precomputeRefs()
{
    std::uint32_t passes = params->getTimeCost();
    std::uint32_t lanes = params->getLanes();
    std::uint32_t segmentBlocks = params->getSegmentBlocks();
    std::uint32_t segmentAddrBlocks =
            (segmentBlocks + ARGON2_REFS_PER_BLOCK - 1)
            / ARGON2_REFS_PER_BLOCK;
    std::uint32_t segments = programContext->getArgon2Type() == ARGON2_ID
            ? lanes * (ARGON2_SYNC_POINTS / 2)
            : passes * lanes * ARGON2_SYNC_POINTS;

    std::size_t shmemSize = THREADS_PER_LANE * sizeof(cl_uint) * 2;

    cl::Kernel kernel = cl::Kernel(programContext->getProgram(),
                                   "argon2_precompute_kernel");
    kernel.setArg<cl::LocalSpaceArg>(0, { shmemSize });
    kernel.setArg<cl::Buffer>(1, refsBuffer);
    kernel.setArg<cl_uint>(2, passes);
    kernel.setArg<cl_uint>(3, lanes);
    kernel.setArg<cl_uint>(4, segmentBlocks);

    cl::NDRange globalRange { THREADS_PER_LANE * segments * segmentAddrBlocks };
    cl::NDRange localRange { THREADS_PER_LANE };
    queue.enqueueNDRangeKernel(kernel, cl::NullRange, globalRange, localRange);
    queue.finish();
}

void KernelRunner::mapMemory(CoinAlgo algo)
{
	if(algo == None)
	    seedHost = queue.enqueueMapBuffer(seedBuffer, true, CL_MAP_WRITE,
                                      0, batchSize * ARGON2_PREHASH_DIGEST_LENGTH);
	else
		seedHost = queue.enqueueMapBuffer(seedBuffer, true, CL_MAP_WRITE,
										  0, 80);

    outHost = queue.enqueueMapBuffer(outBuffer, true, CL_MAP_READ,
                                      0, batchSize * outLen);
}

void KernelRunner::unmapMemory()
{
    queue.enqueueUnmapMemObject(seedBuffer, seedHost);
    queue.enqueueUnmapMemObject(outBuffer, outHost);
}

void KernelRunner::run(CoinAlgo algo, std::uint32_t lanesPerBlock, std::uint32_t jobsPerBlock)
{
    timer = get_time();
    std::uint32_t lanes = params->getLanes();
    std::uint32_t passes = params->getTimeCost();

    if (bySegment) {
        if (lanesPerBlock > lanes || lanes % lanesPerBlock != 0) {
            throw std::logic_error("Invalid lanesPerBlock!");
        }
    } else {
        if (lanesPerBlock != lanes) {
            throw std::logic_error("Invalid lanesPerBlock!");
        }
    }

    if (jobsPerBlock > batchSize || batchSize % jobsPerBlock != 0) {
        throw std::logic_error("Invalid jobsPerBlock!");
    }

    std::size_t shmemSizePreseed = lanes * 2 * 480;

    kpreseed.setArg<cl_uint>(0, algo);
    kpreseed.setArg<cl::LocalSpaceArg>(9, { shmemSizePreseed });
    queue.enqueueNDRangeKernel(kpreseed, cl::NullRange,
                               cl::NDRange(8 * lanes, batchSize), cl::NDRange(8 * lanes, 1));

    cl::NDRange globalRange { THREADS_PER_LANE * lanes, batchSize };
    cl::NDRange localRange { THREADS_PER_LANE * lanesPerBlock, jobsPerBlock };

    std::size_t shmemSize = THREADS_PER_LANE * lanesPerBlock * jobsPerBlock
            * sizeof(cl_uint) * 2;
    kworker.setArg<cl::LocalSpaceArg>(0, { shmemSize });

    if (bySegment) {
        for (std::uint32_t pass = 0; pass < passes; pass++) {
            for (std::uint32_t slice = 0; slice < ARGON2_SYNC_POINTS; slice++) {
                kworker.setArg<cl_uint>(precompute ? 6 : 5, pass);
                kworker.setArg<cl_uint>(precompute ? 7 : 6, slice);
                queue.enqueueNDRangeKernel(kworker, cl::NullRange,
                                           globalRange, localRange);
            }
        }
    } else {
        queue.enqueueNDRangeKernel(kworker, cl::NullRange,
                                   globalRange, localRange);
    }

    std::size_t shmemSizeFinalize = 480;

    kfinalize.setArg<cl::LocalSpaceArg>(5, { shmemSizeFinalize });
    queue.enqueueNDRangeKernel(kfinalize, cl::NullRange,
                               cl::NDRange(THREADS_PER_LANE, batchSize), cl::NDRange(THREADS_PER_LANE, 1));
}

std::uint64_t KernelRunner::finish()
{
    queue.finish();
    return get_time() - timer;
}

uint64_t KernelRunner::get_time() {
#ifdef _WIN64
	FILETIME ft;
	uint64_t tmpres = 0;
	static int tzflag;

	GetSystemTimeAsFileTime(&ft);

	tmpres |= ft.dwHighDateTime;
	tmpres <<= 32;
	tmpres |= ft.dwLowDateTime;

	return tmpres;
#else
	timespec t;
	clock_gettime(CLOCK_MONOTONIC, &t);
	return t.tv_sec * 1000000000 + t.tv_nsec;
#endif
}

void *KernelRunner::getSeedBuffer(int index) {
    return &((uint8_t*)seedHost)[index * ARGON2_PREHASH_DIGEST_LENGTH];
}

void *KernelRunner::getOutBuffer(int index) {
    return &((uint8_t*)outHost)[index * outLen];
}

} // namespace opencl
} // namespace argon2
