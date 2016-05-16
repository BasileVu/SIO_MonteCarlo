#ifndef IMPORTANCE_SAMPLING_H
#define IMPORTANCE_SAMPLING_H

#include "MonteCarloMethod.h"
#include "../generators/RandomValueGenerator.h"

class ImportanceSampling : public MonteCarloMethod {
private:
    InverseFunctions generator;

public:
    ImportanceSampling(const Func& g, const std::vector<double>& xs, const std::vector<double>& ys);

    Sampling sampleWithSize(size_t N);

    Sampling sampleWithMaxWidth(double maxWidth, size_t step);

    Sampling sampleWithMinTime(double maxTime, size_t step);

    void setSeed(const std::seed_seq& seed);

private:
    void sample(size_t step);
};

#endif // IMPORTANCE_SAMPLING_H
