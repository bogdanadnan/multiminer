#include "benchmark.h"

#include <iostream>

int BenchmarkDirector::runBenchmark(Argon2Runner &runner) const
{
    DummyPasswordGenerator pwGen;
    RunTimeStats stats(batchSize);
    for (std::size_t i = 0; i < samples; i++) {
        if (beVerbose) {
            std::cout << "  Sample " << i << "..." << std::endl;
        }
        stats.addSample(runner.runBenchmark(*this, pwGen));
    }
    stats.close();

    if (beVerbose) {
        auto &time = stats.getNanoseconds();
        std::cout << "Mean computation time: "
                  << RunTimeStats::repr((nanosecs)time.getMean())
                  << std::endl;
        std::cout << "Mean deviation: "
                  << RunTimeStats::repr((nanosecs)time.getMeanDeviation())
                  << " (" << time.getMeanDeviationPerMean() * 100.0 << "%)"
                  << std::endl;

        auto &perHash = stats.getNanosecsPerHash();
        std::cout << "Mean computation time (per hash): "
                  << RunTimeStats::repr((nanosecs)perHash.getMean())
                  << std::endl;
        std::cout << "Mean deviation (per hash): "
                  << RunTimeStats::repr((nanosecs)perHash.getMeanDeviation())
                  << std::endl;
        return 0;
    }

    const DataSet *dataSet;
    if (outputType == "ns") {
        dataSet = &stats.getNanoseconds();
    } else if (outputType == "ns-per-hash") {
        dataSet = &stats.getNanosecsPerHash();
    } else {
        std::cerr << progname << ": invalid output type: '"
                  << outputType << "'" << std::endl;
        return 1;
    }

    if (outputMode == "raw") {
        for (auto sample : dataSet->getSamples()) {
            std::cout << sample << std::endl;
        }
    } else if (outputMode == "mean") {
        std::cout << dataSet->getMean() << std::endl;
    } else if (outputMode == "mean-and-mdev") {
        std::cout << dataSet->getMean() << std::endl;
        std::cout << dataSet->getMeanDeviation() << std::endl;
    } else {
        std::cerr << progname << ": invalid output mode: '"
                  << outputMode << "'" << std::endl;
        return 1;
    }
    return 0;
}
