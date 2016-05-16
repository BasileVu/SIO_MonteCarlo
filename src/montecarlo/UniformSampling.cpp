#include <ctime>

#include <vector>
#include "UniformSampling.h"

UniformSampling::UniformSampling(const MonteCarloMethod::Func& g, double a, double b)
        :
        MonteCarloMethod(g), a(a), b(b),
        distribution(std::uniform_real_distribution<double>(0, 1)) {}

MonteCarloMethod::Sampling UniformSampling::sampleWithSize(size_t N) {
    init();

    sample(N);
    return createSampling((double)(clock() - start) / CLOCKS_PER_SEC);
}

MonteCarloMethod::Sampling UniformSampling::sampleWithMaxDelta(double maxDelta, size_t step) {
    init();

    // genere des valeurs tant que la largeur de l'intervalle de confiance est plus grande que "maxDelta"
    do {
        sample(step);
    } while (halfDelta * 2 > maxDelta);

    return createSampling((double)(clock() - start) / CLOCKS_PER_SEC);
}

MonteCarloMethod::Sampling UniformSampling::sampleWithMinTime(double maxTime, size_t step) {
    init();
    double curTime = 0;

    // genere des valeurs tant que le temps maximal d'execution n'est pas atteint
    do {
        clock_t beg = clock();
        sample(step);
        curTime += (double)(clock() - beg) / CLOCKS_PER_SEC;
    } while (curTime < maxTime);

    return createSampling(curTime);
}

void UniformSampling::setSeed(const std::seed_seq &seed) {
    std::seed_seq copy = seed;
    generator.seed(copy);
}

void UniformSampling::sample(size_t step) {
    for (size_t i = 0; i < step; ++i) {
        double X = distribution(generator) * (b - a) + a;
        double Y = g(X);

        sum += Y;
        sumSquares += Y * Y;
    }

    numGen += step;

    mean = sum / numGen;
    double var = (sumSquares / numGen) - mean * mean;
    stdDev = (b-a) * sqrt(var / numGen);
    halfDelta = 1.96 * (b - a) * sqrt(var / numGen);
}

MonteCarloMethod::Sampling UniformSampling::createSampling(double timeElapsed) const {
    double areaEstimator = (b-a) * mean;

    // retourne l'estimateur de l'aire ainsi que l'intervalle de confiance associee (et la taille N, par cohérence)
    return {areaEstimator, stdDev, ConfidenceInterval(areaEstimator, halfDelta), numGen, timeElapsed};
}