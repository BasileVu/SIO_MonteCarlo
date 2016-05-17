#include <ctime>

#include "ImportanceSampling.h"

ImportanceSampling::ImportanceSampling(const std::function<double(double)>& g, const std::vector<double>& xs, const std::vector<double>& ys)
        : MonteCarloMethod(g), generator(xs, ys) {}

MonteCarloMethod::Sampling ImportanceSampling::sampleWithSize(uint64_t N) {
    init();
    sample(N);
    return {mean, stdDev, ConfidenceInterval(mean, halfWidth), N, (double)(clock() - start) / CLOCKS_PER_SEC};
}

MonteCarloMethod::Sampling ImportanceSampling::sampleWithMaxWidth(double maxWidth, uint64_t step) {
    init();

    // genere des valeurs tant que la largeur de l'intervalle de confiance est plus grande que "maxWidth"
    do {
        sample(step);
    } while (halfWidth * 2 > maxWidth);

    return {mean, stdDev, ConfidenceInterval(mean, halfWidth), numGen, (double)(clock() - start) / CLOCKS_PER_SEC};
}

MonteCarloMethod::Sampling ImportanceSampling::sampleWithMinTime(double maxTime, uint64_t step) {
    init();

    double curTime = 0;

    // genere des valeurs tant que le temps maximal d'execution n'est pas atteint
    do {
        clock_t beg = clock();
        sample(step);
        curTime += (double)(clock() - beg) / CLOCKS_PER_SEC;
    } while (curTime < maxTime);

    return {mean, stdDev, ConfidenceInterval(mean, halfWidth), numGen, curTime};
}

void ImportanceSampling::setSeed(const std::seed_seq &seed) {
    generator.setSeed(seed);
}

void ImportanceSampling::sample(uint64_t step) {
    const PiecewiseLinearFunction& f = generator.getPWLFunc();

    for (uint64_t i = 0; i < step; ++i) {
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
    halfWidth = 1.96 * stdDev;
}


