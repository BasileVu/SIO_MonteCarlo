#include <ctime>

#include <vector>
#include "UniformSampling.h"

UniformSampling::UniformSampling(const MonteCarloMethod::Func& g, double a, double b)
        :
        MonteCarloMethod(g), a(a), b(b),
        distribution(std::uniform_real_distribution<double>(0, 1)) {}

MonteCarloMethod::Sampling UniformSampling::sampleWithSize(size_t N) {

    clock_t start = clock();

    // genere N valeurs et retourne dans res la moyenne et la demi-largeur de l'IC
    double S = 0, Q = 0;
    size_t tmpN = 0; // sauvegarde de N

    Result res = sample(N, tmpN, S, Q);
    return createSampling(res.mean, res.halfDelta, N, (double)(clock() - start) / CLOCKS_PER_SEC);
}

MonteCarloMethod::Sampling UniformSampling::sampleWithMaxDelta(double maxDelta, size_t step) {

    clock_t start = clock();

    double S = 0, Q = 0;
    size_t N = 0;

    // genere des valeurs tant que la largeur de l'intervalle de confiance est plus grande que "maxDelta"
    Result res;
    do {
        res = sample(step, N, S, Q);
    } while (res.halfDelta * 2 > maxDelta);

    return createSampling(res.mean, res.halfDelta, N, (double)(clock() - start) / CLOCKS_PER_SEC);
}

MonteCarloMethod::Sampling UniformSampling::sampleWithMaxTime(double maxTime, size_t step) {

    double S = 0, Q = 0;
    size_t N = 0;
    double curTime = 0;

    // genere des valeurs tant que le temps maximal d'execution n'est pas atteint
    Result res;
    do {
        clock_t beg = clock();
        res = sample(step, N, S, Q);
        curTime += (double)(clock() - beg) / CLOCKS_PER_SEC;
    } while (curTime < maxTime);

    return createSampling(res.mean, res.halfDelta, N, curTime);
}

void UniformSampling::setSeed(const std::seed_seq &seed) {
    std::seed_seq copy = seed;
    generator.seed(copy);
}


UniformSampling::Result UniformSampling::sample(size_t step, size_t& N, double& S, double& Q) {

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

MonteCarloMethod::Sampling UniformSampling::createSampling(double mean, double halfDelta, size_t N, double timeElapsed) const {
    double areaEstimator = (b-a) * mean;

    // retourne l'estimateur de l'aire ainsi que l'intervalle de confiance associee (et la taille N, par cohérence)
    return {areaEstimator, ConfidenceInterval(areaEstimator, halfDelta), N, timeElapsed};
}