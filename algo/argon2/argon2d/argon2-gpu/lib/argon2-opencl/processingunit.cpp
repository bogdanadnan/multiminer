#include "processingunit.h"

#include <limits>
#ifndef NDEBUG
#include <iostream>
#endif

namespace argon2 {
namespace opencl {

static bool isPowerOfTwo(std::uint32_t x)
{
    return (x & (x - 1)) == 0;
}

ProcessingUnit::ProcessingUnit(
        const ProgramContext *programContext, const Argon2Params *params,
        const Device *device, std::size_t batchSize, const CoinAlgo algo,
		bool bySegment, bool precomputeRefs)
    : programContext(programContext), params(params), device(device),
      runner(programContext, params, device, batchSize, params->getOutputLength(), bySegment,
             precomputeRefs, (std::uint8_t*)params->getSecret(), params->getSecretLength(),
			 (std::uint8_t*)params->getAssocData(), params->getAssocDataLength()),
      bestLanesPerBlock(runner.getMinLanesPerBlock()),
      bestJobsPerBlock(runner.getMinJobsPerBlock()),
      algo(algo)
{
    runner.mapMemory(None);
    /* pre-fill first blocks with pseudo-random data: */
    for (std::size_t i = 0; i < batchSize; i++) {
        setPassword(i, NULL, 0);
    }
    runner.unmapMemory();

    if (runner.getMaxLanesPerBlock() > runner.getMinLanesPerBlock()
            && isPowerOfTwo(runner.getMaxLanesPerBlock())) {
#ifndef NDEBUG
        std::cerr << "[INFO] Tuning lanes per block..." << std::endl;
#endif

        std::uint64_t bestTime = 0xFFFFFFFFFFFFFFFF;
        for (std::uint32_t lpb = 1; lpb <= runner.getMaxLanesPerBlock();
             lpb *= 2)
        {
            std::uint64_t time;
            try {
                runner.run(None, lpb, bestJobsPerBlock);
                time = runner.finish();
            } catch(cl::Error &ex) {
#ifndef NDEBUG
                std::cerr << "[WARN]   OpenCL error on " << lpb
                          << " lanes per block: " << ex.what() << std::endl;
#endif
                break;
            }

#ifndef NDEBUG
            std::cerr << "[INFO]   " << lpb << " lanes per block: "
                      << time/1000000.0 << " ms" << std::endl;
#endif

            if (time < bestTime) {
                bestTime = time;
                bestLanesPerBlock = lpb;
            }
        }
#ifndef NDEBUG
        std::cerr << "[INFO] Picked " << bestLanesPerBlock
                  << " lanes per block." << std::endl;
#endif
    }

    /* Only tune jobs per block if we hit maximum lanes per block: */
    if (bestLanesPerBlock == runner.getMaxLanesPerBlock()
            && runner.getMaxJobsPerBlock() > runner.getMinJobsPerBlock()
            && isPowerOfTwo(runner.getMaxJobsPerBlock())) {
#ifndef NDEBUG
        std::cerr << "[INFO] Tuning jobs per block..." << std::endl;
#endif

		std::uint64_t bestTime = 0xFFFFFFFFFFFFFFFF;
        for (std::uint32_t jpb = 1; jpb <= runner.getMaxJobsPerBlock();
             jpb *= 2)
        {
			std::uint64_t time;
            try {
                runner.run(None, bestLanesPerBlock, jpb);
                time = runner.finish();
            } catch(cl::Error &ex) {
#ifndef NDEBUG
                std::cerr << "[WARN]   OpenCL error on " << jpb
                          << " jobs per block: " << ex.what() << std::endl;
#endif
                break;
            }

#ifndef NDEBUG
            std::cerr << "[INFO]   " << jpb << " jobs per block: "
                      << time/1000000.0 << " ms" << std::endl;
#endif

            if (time < bestTime) {
                bestTime = time;
                bestJobsPerBlock = jpb;
            }
        }
#ifndef NDEBUG
        std::cerr << "[INFO] Picked " << bestJobsPerBlock
                  << " jobs per block." << std::endl;
#endif
    }

    runner.mapMemory(algo);
}

void ProcessingUnit::setPassword(std::size_t index, const void *pw,
                                 std::size_t pwSize)
{
    void *memory = runner.getSeedBuffer(index);
    params->fillFirstBlocks(memory, pw, pwSize,
                            programContext->getArgon2Type(),
                            programContext->getArgon2Version());
}

void ProcessingUnit::setPasswordSameSalt(std::size_t index, const void *pw,
                                 std::size_t pwSize)
{
    void *memory = runner.getSeedBuffer(index);
    params->fillFirstBlocksSameSalt(memory, pw, pwSize,
                            programContext->getArgon2Type(),
                            programContext->getArgon2Version());
}

void ProcessingUnit::setInput(const void *pw,
                              std::size_t pwSize)
{
    void *buffer = runner.getSeedBuffer(0);
    memcpy(buffer, pw, pwSize);
}

void *ProcessingUnit::getHash(std::size_t index)
{
	return runner.getOutBuffer(index);
}

void ProcessingUnit::beginProcessing()
{
	runner.unmapMemory();
    runner.run(algo, bestLanesPerBlock, bestJobsPerBlock);
}

void ProcessingUnit::endProcessing()
{
	runner.mapMemory(algo);
    runner.finish();
}

} // namespace opencl
} // namespace argon2
