#ifndef UNIFORM_SAMPLING_H
#define UNIFORM_SAMPLING_H

#include "MonteCarloMethod.h"

class UniformSampling : public MonteCarloMethod {
private:
    std::mt19937_64 generator;
    std::uniform_real_distribution<double> distribution;

public:
    UniformSampling(const Func& g);

    Sampling sample(size_t numPoints, double a, double b);

    void setSeed(const std::seed_seq& seed);
};

#endif // UNIFORMS_AMPLING_H
