#ifndef UNIFORM_SAMPLING_H
#define UNIFORM_SAMPLING_H

#include "MonteCarloMethod.h"

class UniformSampling : public MonteCarloMethod {
private:
    std::mt19937_64 mtGenerator;
    std::uniform_real_distribution<double> uniformDistr;

    double a, b;

public:
    UniformSampling(const Func& g, double a, double b);

    Sampling sampleWithSize(size_t N);

    Sampling sampleWithMaxWidth(double maxWidth, size_t step);

    Sampling sampleWithMinTime(double maxTime, size_t step);

    void setSeed(const std::seed_seq& seed);

private:
    void sample(size_t step);

    Sampling createSampling(double timeElapsed) const;
};

#endif // UNIFORMS_AMPLING_H
