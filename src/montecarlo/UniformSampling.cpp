#include <vector>
#include "UniformSampling.h"

UniformSampling::UniformSampling(const MonteCarloMethod::Func& g)
        : MonteCarloMethod(g), distribution(std::uniform_real_distribution<double>(0, 1)) {}

MonteCarloMethod::Sampling UniformSampling::sample(size_t N, double a, double b) {
    double sum = 0, sumSquares = 0;

    for (size_t i = 0; i < N; ++i) {
        double X = distribution(generator) * (b-a) + a;
        double Y = g(X);

        sum += Y;
        sumSquares += Y*Y;
    }

    double mean = sum / N;
    double var = sumSquares / N - mean*mean;
    double areaEstimator = (b-a) * mean;

    double halfDelta = 1.96 * (b-a) * sqrt(var/N);

    // on retourne l'estimateur de l'aire ainsi que l'intervalle de confiance associee.
    return {areaEstimator, {areaEstimator - halfDelta, areaEstimator + halfDelta, halfDelta * 2}};
}

void UniformSampling::setSeed(const std::seed_seq &seed) {
    std::seed_seq copy = seed;
    generator.seed(copy);
}

