#ifndef CONTROLVARIABLE_H
#define CONTROLVARIABLE_H

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
    struct Result {
        double SV;
        double QV;
        double c;
    };

    // TODO
    Result firstStep(size_t M, double a, double b, const PiecewiseLinearFunction &h, double mu);
};

#endif // CONTROLVARIABLE_H
