#include <vector>
#include "UniformSampling.h"

UniformSampling::UniformSampling(const MonteCarloMethod::Func& g)
        : MonteCarloMethod(g), distribution(std::uniform_real_distribution<double>(0, 1)) {}

MonteCarloMethod::Sampling UniformSampling::sample(size_t N, double a, double b) {

    // genere N valeurs et retourne dans res la moyenne et la demi-largeur de l'IC
    double S = 0, Q = 0;
    Result res = sampleN(N, N, a, b, S, Q);

    double areaEstimator = (b-a) * res.mean;

    // retourne l'estimateur de l'aire ainsi que l'intervalle de confiance associee (et la taille N, par cohérence)
    return {areaEstimator, ConfidenceInterval(areaEstimator, res.halfDelta), N};
}

MonteCarloMethod::Sampling UniformSampling::sample(double maxDelta, size_t step, double a, double b) {

    double S = 0, Q = 0;
    size_t N = 0;

    // genere des valeurs tant que la largeur de l'intervalle de confiance est plus grande que "maxDelta"
    Result res;
    do {
        res = sampleN(step, N, a, b, S, Q);
    } while (res.halfDelta * 2 > maxDelta);

    double areaEstimator = (b - a) * res.mean;

    // retourne l'estimateur de l'aire, l'intervalle de confiance associee ainsi que la taille N
    return {areaEstimator, ConfidenceInterval(areaEstimator, res.halfDelta), N};
}

UniformSampling::Result UniformSampling::sampleN(size_t step, size_t& N, double a, double b, double& S, double& Q) {
    for (size_t i = 0; i < step; ++i) {
        double X = distribution(generator) * (b - a) + a;
        double Y = g(X);

        S += Y;
        Q += Y * Y;
    }

    N += step;

    double mean = S / N;
    double var = Q / N - mean * mean;
    double halfDelta = 1.96 * (b - a) * sqrt(var / N);

    return {mean, halfDelta};
}

void UniformSampling::setSeed(const std::seed_seq &seed) {
    std::seed_seq copy = seed;
    generator.seed(copy);
}