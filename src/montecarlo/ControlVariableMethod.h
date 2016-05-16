#ifndef CONTROL_VARIABLE_H
#define CONTROL_VARIABLE_H

#include "MonteCarloMethod.h"

class ControlVariable : public MonteCarloMethod {
private:
    std::mt19937_64 generator;
    std::uniform_real_distribution<double> distribution;

    PiecewiseLinearFunction h;

    double a, b;
    double mu;

public:
    ControlVariable(const Func& g, double a, double b, const std::vector<double>& xs, const std::vector<double>& ys);

    Sampling sampleWithSize(size_t M, size_t N);

    Sampling sampleWithMaxDelta(size_t M, double maxDelta, size_t step);

    Sampling sampleWithMinTime(size_t M, double maxTime, size_t step);

    void setSeed(const std::seed_seq& seed);

private:
    struct ResultFirst {
        double SV;
        double QV;
        double c;
    };

    // TODO
    ResultFirst firstStep(size_t M);

    void secondStep(size_t step, size_t& N, double c, double& SV, double& QV);

    Sampling createSampling(size_t N, double timeElapsed);
};

#endif // CONTROL_VARIABLE_H