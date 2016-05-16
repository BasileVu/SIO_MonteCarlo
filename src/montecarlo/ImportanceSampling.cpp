#include <ctime>

#include "ImportanceSampling.h"

ImportanceSampling::ImportanceSampling(const std::function<double(double)>& g, const std::vector<double>& xs, const std::vector<double>& ys)
        : MonteCarloMethod(g), generator(xs, ys) {}

MonteCarloMethod::Sampling ImportanceSampling::sampleWithSize(size_t N) {
    init();
    sample(N);
    return {mean, stdDev, ConfidenceInterval(mean, halfDelta), N, (double)(clock() - start) / CLOCKS_PER_SEC};
}

MonteCarloMethod::Sampling ImportanceSampling::sampleWithMaxDelta(double maxDelta, size_t step) {
    init();

    // genere des valeurs tant que la largeur de l'intervalle de confiance est plus grande que "maxDelta"
    do {
        sample(step);
    } while (halfDelta * 2 > maxDelta);

    return {mean, stdDev, ConfidenceInterval(mean, halfDelta), numGen, (double)(clock() - start) / CLOCKS_PER_SEC};
}

MonteCarloMethod::Sampling ImportanceSampling::sampleWithMinTime(double maxTime, size_t step) {
    init();

    double curTime = 0;

    // genere des valeurs tant que le temps maximal d'execution n'est pas atteint
    do {
        clock_t beg = clock();
        sample(step);
        curTime += (double)(clock() - beg) / CLOCKS_PER_SEC;
    } while (curTime < maxTime);

    return {mean, stdDev, ConfidenceInterval(mean, halfDelta), numGen, curTime};
}

void ImportanceSampling::setSeed(const std::seed_seq &seed) {
    generator.setSeed(seed);
}

void ImportanceSampling::sample(size_t step) {
    const PiecewiseLinearFunction& f = generator.getPWLFunc();

    for (size_t i = 0; i < step; ++i) {
        double X = generator.generate();
        double Y = g(X) / f(X);

        sum += Y;
        sumSquares += Y*Y;
    }

    numGen += step;

    // multiplication a la fin plutot que multiplier Y a chaque iteration dans la boucle
    double tmpS = sum * f.A;
    double tmpQ = sumSquares * (f.A * f.A);

    mean = tmpS/numGen;
    double var = tmpQ/numGen - mean*mean;
    stdDev = sqrt(var/numGen);
    halfDelta = 1.96 * stdDev;
}


