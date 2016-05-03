#ifndef CONTROL_VARIABLE_H
#define CONTROL_VARIABLE_H

#include "MonteCarloMethod.h"

class ControlVariable : public MonteCarloMethod {
private:
    std::mt19937_64 generator;
    std::uniform_real_distribution<double> distribution;

public:
    ControlVariable(const Func &g);

    Sampling sample(size_t M, size_t N, double a, double b,
                    const std::vector<double>& xs, const std::vector<double>& ys);

    Sampling sample(size_t M, double maxDelta, size_t step, double a, double b,
                    const std::vector<double>& xs, const std::vector<double>& ys);

    void setSeed(const std::seed_seq& seed);

private:
    struct ResultFirst {
        double SV;
        double QV;
        double c;
    };

    struct ResultSecond {
        double meanV;
        double halfDelta;
    };

    // TODO
    ResultFirst firstStep(size_t M, double a, double b, const PiecewiseLinearFunction &h, double mu);

    ResultSecond secondStep(size_t step, size_t& N, double a, double b, const PiecewiseLinearFunction& h,
                                            double mu, double c, double& SV, double& QV);
};

#endif // CONTROL_VARIABLE_H
