#ifndef MONTECARLOMETHOD_H
#define MONTECARLOMETHOD_H

#include <random>
#include <functional>
#include "Stats.h"

class MonteCarloMethod {
protected:
    const std::function<double(double)>& func;
    std::mt19937_64 generator;
    std::uniform_real_distribution<double> distribution;

public:

    struct Sampling {
        double areaEstimator;
        ConfidenceInterval confidenceInterval;
    };

    MonteCarloMethod(const std::function<double(double)>& func);
    virtual Sampling sample() = 0;
};

#endif // MONTECARLOMETHOD_H
