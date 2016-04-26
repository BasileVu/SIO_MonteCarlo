#include "ImportanceSampling.h"
#include "RandomValueGenerator.h"

ImportanceSampling::ImportanceSampling(const std::function<double(double)>& g)
        : MonteCarloMethod(g) {
}

MonteCarloMethod::Sampling ImportanceSampling::sample(size_t numPoints, const MonteCarloMethod::Func& f,
                                                      const std::vector<double> xs, const std::vector<double> ys) {

    double sum = 0, squares = 0;

    InverseFunctions inv(xs, ys);
    inv.setSeed(seed);

    for (size_t i = 0; i < numPoints; ++i) {
        double X = inv.generate();
        double Y = g(X) / f(X);

        sum += Y;
        squares += Y*Y;
    }

    double mean = sum / numPoints;
    double var = squares/numPoints - mean*mean;

    double halfDelta = 1.96 * sqrt(var/numPoints);

    return {mean, {mean - halfDelta, mean + halfDelta}};
}

void ImportanceSampling::setSeed(const std::seed_seq &seed) {
    this->seed = seed;
}

