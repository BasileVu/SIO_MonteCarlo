#ifndef IMPORTANCE_SAMPLING_H
#define IMPORTANCE_SAMPLING_H

#include "MonteCarloMethod.h"
#include "../generators/RandomValueGenerator.h"

class ImportanceSampling : public MonteCarloMethod {
private:
    InverseFunctions generator;

public:
    ImportanceSampling(const Func& g, const std::vector<double>& xs, const std::vector<double>& ys);

    Sampling sampleWithSize(size_t N);

    Sampling sampleWithMaxDelta(double maxDelta, size_t step);

    Sampling sampleWithMaxTime(double maxTime, size_t step);

    void setSeed(const std::seed_seq& seed);

private:
    struct Result {
        double mean;
        double halfDelta;
    };

    Result sample(size_t step, size_t& N, double& S, double& Q);
};

#endif // IMPORTANCE_SAMPLING_H
