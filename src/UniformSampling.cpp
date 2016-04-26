#include <vector>
#include "UniformSampling.h"

UniformSampling::UniformSampling(const MonteCarloMethod::Func& g)
        : MonteCarloMethod(g), distribution(std::uniform_real_distribution<double>(0, 1)) {}

MonteCarloMethod::Sampling UniformSampling::sample(size_t numPoints) {

    const double a = 0, b = 15;
    double sum = 0, sumSquares = 0;

    for (size_t i = 0; i < numPoints; ++i) {
        double X = distribution(generator) * (b-a);
        double Y = g(X);

        sum += Y;
        sumSquares += Y*Y;
    }

    double mean = sum / numPoints;
    double var = sumSquares / numPoints - mean*mean;
    double areaEstimator = (b-a) * mean;

    double halfDelta = 1.96 * (b-a) * sqrt(var/numPoints);

    // on retourne l'estimateur de l'aire ainsi que l'intervalle de confiance associee.
    return {areaEstimator, {areaEstimator - halfDelta, areaEstimator + halfDelta, halfDelta * 2}};
}

void UniformSampling::setSeed(const std::seed_seq &seed) {
    std::seed_seq copy = seed;
    generator.seed(copy);
}

