#ifndef IMPORTANCE_SAMPLING_H
#define IMPORTANCE_SAMPLING_H

#include "MonteCarloMethod.h"
#include "../generators/RandomValueGenerator.h"

class ImportanceSampling : public MonteCarloMethod {
private:
    std::seed_seq seed;

public:
    ImportanceSampling(const Func& g);

    Sampling sample(double maxDelta, size_t step, const std::vector<double> xs, const std::vector<double> ys);

    Sampling sample(size_t N, const std::vector<double> xs, const std::vector<double> ys);

    void setSeed(const std::seed_seq& seed);

private:
    struct Result {
        double mean;
        double halfDelta;
    };

    Result sample(size_t step, size_t& N, RandomValueGenerator& rvg, double& S, double& Q);
};

#endif // IMPORTANCE_SAMPLING_H
