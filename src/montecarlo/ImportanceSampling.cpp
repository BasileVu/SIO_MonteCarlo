#include "ImportanceSampling.h"
#include "../generators/RandomValueGenerator.h"

ImportanceSampling::ImportanceSampling(const std::function<double(double)>& g)
        : MonteCarloMethod(g) {
}

MonteCarloMethod::Sampling ImportanceSampling::sample(size_t N, const std::vector<double> xs, const std::vector<double> ys) {

    double sum = 0, squares = 0;

    InverseFunctions inv(xs, ys);
    inv.setSeed(seed);
    const PiecewiseLinearFunction& f = inv.getPWLFunc();

    for (size_t i = 0; i < N; ++i) {
        double X = inv.generate();
        double Y = g(X) / (f.pieces[f.findPart(X)].f_k(X) / f.A);

        sum += Y;
        squares += Y*Y;
    }

    double mean = sum / N;
    double var = squares/N - mean*mean;

    double halfDelta = 1.96 * sqrt(var/N);

    return {mean, {mean - halfDelta, mean + halfDelta}};
}

void ImportanceSampling::setSeed(const std::seed_seq &seed) {
    this->seed = seed;
}

