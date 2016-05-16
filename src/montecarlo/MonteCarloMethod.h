#ifndef MONTECARLOMETHOD_H
#define MONTECARLOMETHOD_H

#include <random>
#include <functional>
#include <ctime>

#include "../utility/Stats.h"

class MonteCarloMethod {
public:
    typedef std::function<double(double)> Func;

protected:
    const Func& g;

    double mean;
    double stdDev;
    double halfDelta;

    double sum;
    double sumSquares;
    size_t numGen;

    clock_t start;

public:
    struct Sampling {
        double areaEstimator;                   // aire estimee
        double stdDevEstimator;                 // estimateur de l'ecart-type de l'aire estimee
        ConfidenceInterval confidenceInterval;  // intervalle de confiance a 95%
        size_t N;                               // taille de l'echantillon
        double elapsedTime;                     // temps pour creer la totalite de l'echantillon
    };

    MonteCarloMethod(const Func& g);
    virtual void setSeed(const std::seed_seq& seed) = 0;

protected:
    void init();
};

#endif // MONTECARLOMETHOD_H
