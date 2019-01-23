#ifndef ARGON2_OPENCL_KERNELRUNNER_H
#define ARGON2_OPENCL_KERNELRUNNER_H

#include "programcontext.h"
#include "argon2-gpu-common/argon2params.h"

namespace argon2 {
namespace opencl {

class KernelRunner
{
private:
    const ProgramContext *programContext;
    const Argon2Params *params;

    std::size_t batchSize, outLen;
    std::uint64_t timer;
    bool bySegment;
    bool precompute;

    cl::CommandQueue queue;
	cl::Kernel kpreseed, kworker, kfinalize;
    cl::Buffer memoryBuffer, refsBuffer, seedBuffer, outBuffer, secretBuffer, adBuffer;
    void *seedHost, *outHost;
    std::size_t secretLen, adLen;

    std::size_t memorySize;

    void precomputeRefs();
	uint64_t get_time();

public:
    std::uint32_t getMinLanesPerBlock() const
    {
        return bySegment ? 1 : params->getLanes();
    }
    std::uint32_t getMaxLanesPerBlock() const { return params->getLanes(); }

    std::uint32_t getMinJobsPerBlock() const { return 1; }
    std::uint32_t getMaxJobsPerBlock() const { return batchSize; }

    std::size_t getBatchSize() const { return batchSize; }
    void *getSeedBuffer(int index);
    void *getOutBuffer(int index);

    KernelRunner(const ProgramContext *programContext,
                 const Argon2Params *params, const Device *device,
                 std::size_t batchSize, std::size_t outLen,
                 bool bySegment, bool precompute,
                 std::uint8_t *secret, std::size_t secretLen,
				 std::uint8_t *ad, std::size_t adLen);

	void mapMemory(CoinAlgo algo);
	void unmapMemory();

    void run(CoinAlgo algo, std::uint32_t lanesPerBlock, std::uint32_t jobsPerBlock);
    uint64_t finish();
};

} // namespace opencl
} // namespace argon2

#endif // ARGON2_OPENCL_KERNELRUNNER_H
