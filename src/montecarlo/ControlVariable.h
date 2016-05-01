#ifndef CONTROLVARIABLE_H
#define CONTROLVARIABLE_H

#include "MonteCarloMethod.h"

class ControlVariable : public MonteCarloMethod {
private:
    std::mt19937_64 generator;
    std::uniform_real_distribution<double> distribution;

public:
    ControlVariable(const Func &g);

    Sampling sample(size_t M, size_t N, double a, double b, const Func& h, double mu);

    void setSeed(const std::seed_seq& seed);
};

#endif // CONTROLVARIABLE_H
