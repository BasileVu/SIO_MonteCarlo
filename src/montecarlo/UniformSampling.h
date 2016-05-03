#ifndef UNIFORM_SAMPLING_H
#define UNIFORM_SAMPLING_H

#include "MonteCarloMethod.h"

class UniformSampling : public MonteCarloMethod {
private:
    std::mt19937_64 generator;
    std::uniform_real_distribution<double> distribution;

    double a, b;

public:
    UniformSampling(const Func& g, double a, double b);

    Sampling sample(size_t N);

    Sampling sample(double maxDelta, size_t step);

    void setSeed(const std::seed_seq& seed);

private:
    struct Result {
        double mean;
        double halfDelta;
    };

    Result sample(size_t step, size_t &N, double &S, double &Q);

    Sampling createSampling(double mean, double halfDelta, size_t N) const;
};

#endif // UNIFORMS_AMPLING_H
