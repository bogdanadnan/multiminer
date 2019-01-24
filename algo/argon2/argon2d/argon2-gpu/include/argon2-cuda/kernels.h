#ifndef ARGON2_CUDA_KERNELS_H
#define ARGON2_CUDA_KERNELS_H

#if HAVE_CUDA

#include <cuda_runtime.h>
#include <cstdint>
#include <vector>

/* workaround weird CMake/CUDA bug: */
#ifdef argon2
#undef argon2
#endif

namespace argon2 {
namespace cuda {

class KernelRunner
{
private:
    std::uint32_t type, version;
    std::uint32_t passes, lanes, segmentBlocks;
    std::size_t batchSize, outLen;
    std::uint64_t timer;
    bool bySegment;
    bool precompute;

    cudaStream_t stream;
	void *memory;
	void *seed;
	void *seed_host;
	void *out;
	void *out_host;
    void *refs;
    std::uint32_t *secret;
	std::size_t secretLen;
    std::uint32_t *ad;
    std::size_t adLen;

    void precomputeRefs();

    void runKernelSegment(std::uint32_t lanesPerBlock,
                          std::uint32_t jobsPerBlock,
                          std::uint32_t pass, std::uint32_t slice);
    void runKernelOneshot(std::uint32_t lanesPerBlock,
                          std::uint32_t jobsPerBlock);

	void runKernelPreseed(CoinAlgo algo);
	void runKernelFinalize();

    void synchronize();

    uint64_t get_time();

public:
    std::uint32_t getMinLanesPerBlock() const { return bySegment ? 1 : lanes; }
    std::uint32_t getMaxLanesPerBlock() const { return lanes; }

    std::uint32_t getMinJobsPerBlock() const { return 1; }
    std::uint32_t getMaxJobsPerBlock() const { return batchSize; }

    std::size_t getBatchSize() const { return batchSize; }
	void *getSeedBuffer(int index);
	void *getOutBuffer(int index);

    KernelRunner(std::uint32_t type, std::uint32_t version,
                 std::uint32_t passes, std::uint32_t lanes,
                 std::uint32_t segmentBlocks, std::size_t batchSize, std::size_t outLen,
                 std::int32_t deviceIndex,
                 bool bySegment, bool precompute,
                 std::uint8_t *secret, std::size_t secretLen,
                 std::uint8_t *ad, std::size_t adLen);
    ~KernelRunner();

    void writeInputMemory(CoinAlgo algo);
    void readOutputMemory();

    void run(CoinAlgo algo, std::uint32_t lanesPerBlock, std::uint32_t jobsPerBlock);
    uint64_t finish();
};

} // cuda
} // argon2

#endif /* HAVE_CUDA */

#endif // ARGON2_CUDA_KERNELS_H
