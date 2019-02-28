#include "processingunit.h"

#include "cudaexception.h"

#include <limits>
#ifndef NDEBUG
#include <iostream>
#endif
#include <cstring>

namespace argon2 {
namespace cuda {

static void setCudaDevice(int deviceIndex)
{
    int currentIndex = -1;
    CudaException::check(cudaGetDevice(&currentIndex));
    if (currentIndex != deviceIndex) {
        CudaException::check(cudaSetDevice(deviceIndex));
    }
}

static bool isPowerOfTwo(std::uint32_t x)
{
    return (x & (x - 1)) == 0;
}

ProcessingUnit::ProcessingUnit(
        const ProgramContext *programContext, const Argon2Params *params,
        const Device *device, std::size_t batchSize, const CoinAlgo algo, bool bySegment,
        bool precomputeRefs)
    : programContext(programContext), params(params), device(device),
      runner(programContext->getArgon2Type(),
             programContext->getArgon2Version(), params->getTimeCost(),
             params->getLanes(), params->getSegmentBlocks(), batchSize, params->getOutputLength(),
             device->getDeviceIndex(),
             bySegment, precomputeRefs,
             (std::uint8_t*)params->getSecret(), params->getSecretLength(),
             (std::uint8_t*)params->getAssocData(), params->getAssocDataLength()),
      bestLanesPerBlock(runner.getMinLanesPerBlock()),
      bestJobsPerBlock(runner.getMinJobsPerBlock()),
      algo(algo)
{
    setCudaDevice(device->getDeviceIndex());

    /* pre-fill first blocks with pseudo-random data: */
    for (std::size_t i = 0; i < batchSize; i++) {
        setPassword(i, NULL, 0);
    }

    runner.writeInputMemory(None);

    if (runner.getMaxLanesPerBlock() > runner.getMinLanesPerBlock()
            && isPowerOfTwo(runner.getMaxLanesPerBlock())) {
#ifndef NDEBUG
        std::cerr << "[INFO] Tuning lanes per block..." << std::endl;
#endif

        uint64_t bestTime = 0xFFFFFFFFFFFFFFFF;
        for (std::uint32_t lpb = 1; lpb <= runner.getMaxLanesPerBlock();
             lpb *= 2)
        {
            uint64_t time;
            try {
                runner.run(None, lpb, bestJobsPerBlock);
                time = runner.finish();
            } catch(CudaException &ex) {
#ifndef NDEBUG
                std::cerr << "[WARN]   CUDA error on " << lpb
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

        uint64_t bestTime = 0xFFFFFFFFFFFFFFFF;
        for (std::uint32_t jpb = 1; jpb <= runner.getMaxJobsPerBlock();
             jpb *= 2)
        {
            uint64_t time;
            try {
                runner.run(None, bestLanesPerBlock, jpb);
                time = runner.finish();
            } catch(CudaException &ex) {
#ifndef NDEBUG
                std::cerr << "[WARN]   CUDA error on " << jpb
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
}

void ProcessingUnit::setPassword(std::size_t index, const void *pw,
                                 std::size_t pwSize)
{
	void *buffer = runner.getSeedBuffer(index);
    params->fillFirstBlocks(buffer, pw, pwSize,
                            programContext->getArgon2Type(),
                            programContext->getArgon2Version());
}

void ProcessingUnit::setPasswordSameSalt(std::size_t index, const void *pw,
                                 std::size_t pwSize)
{
    void *buffer = runner.getSeedBuffer(index);
    params->fillFirstBlocksSameSalt(buffer, pw, pwSize,
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
    setCudaDevice(device->getDeviceIndex());
    runner.writeInputMemory(algo);
    runner.run(algo, bestLanesPerBlock, bestJobsPerBlock);
}

void ProcessingUnit::endProcessing()
{
	runner.readOutputMemory();
	runner.finish();
}

} // namespace cuda
} // namespace argon2
