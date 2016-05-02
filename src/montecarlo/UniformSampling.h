#ifndef UNIFORM_SAMPLING_H
#define UNIFORM_SAMPLING_H

#include "MonteCarloMethod.h"

class UniformSampling : public MonteCarloMethod {
private:
    std::mt19937_64 generator;
    std::uniform_real_distribution<double> distribution;

public:
    UniformSampling(const Func& g);

    Sampling sample(size_t N, double a, double b);

    Sampling sample(double maxDelta, size_t step, double a, double b);

    void setSeed(const std::seed_seq& seed);

private:
    struct Result {
        double mean;
        double halfDelta;
    };

    Result sampleN(size_t step, size_t& N, double a, double b, double& S, double& Q);
};

#endif // UNIFORMS_AMPLING_H
