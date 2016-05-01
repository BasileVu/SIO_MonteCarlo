#ifndef IMPORTANCE_SAMPLING_H
#define IMPORTANCE_SAMPLING_H

#include "MonteCarloMethod.h"

class ImportanceSampling : public MonteCarloMethod {
private:
    std::seed_seq seed;

public:
    ImportanceSampling(const Func& g);

    Sampling sample(size_t N, const std::vector<double> xs, const std::vector<double> ys);

    void setSeed(const std::seed_seq& seed);
};

#endif // IMPORTANCE_SAMPLING_H
