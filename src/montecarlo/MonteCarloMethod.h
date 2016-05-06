#ifndef MONTECARLOMETHOD_H
#define MONTECARLOMETHOD_H

#include <random>
#include <functional>
#include <ctime>

#include "../utility/Stats.h"

class MonteCarloMethod {
public:
    typedef std::function<double(double)>Func;

protected:
    const Func& g;

public:

    struct Sampling {
        double areaEstimator;
        ConfidenceInterval confidenceInterval;
        size_t N;
        double elapsedTime;
    };

    MonteCarloMethod(const Func& g);
    virtual void setSeed(const std::seed_seq& seed) = 0;
};

#endif // MONTECARLOMETHOD_H
