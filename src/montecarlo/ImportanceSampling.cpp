#include <ctime>

#include "ImportanceSampling.h"

ImportanceSampling::ImportanceSampling(const std::function<double(double)>& g, const std::vector<double>& xs, const std::vector<double>& ys)
        : MonteCarloMethod(g), generator(xs, ys) {}

MonteCarloMethod::Sampling ImportanceSampling::sampleWithSize(size_t N) {

    clock_t start = clock();
    double S = 0, Q = 0;
    size_t tmpN = 0;

    sample(N, tmpN, S, Q);

    return {mean, stdDev, ConfidenceInterval(mean, halfDelta), N, (double)(clock() - start) / CLOCKS_PER_SEC};
}

MonteCarloMethod::Sampling ImportanceSampling::sampleWithMaxDelta(double maxDelta, size_t step) {

    clock_t start = clock();
    double S = 0, Q = 0;
    size_t N = 0;

    // genere des valeurs tant que la largeur de l'intervalle de confiance est plus grande que "maxDelta"
    do {
        sample(step, N, S, Q);
    } while (halfDelta * 2 > maxDelta);

    return {mean, stdDev, ConfidenceInterval(mean, halfDelta), N, (double)(clock() - start) / CLOCKS_PER_SEC};
}

MonteCarloMethod::Sampling ImportanceSampling::sampleWithMinTime(double maxTime, size_t step) {

    double S = 0, Q = 0;
    size_t N = 0;
    double curTime = 0;

    // genere des valeurs tant que le temps maximal d'execution n'est pas atteint
    do {
        clock_t beg = clock();
        sample(step, N, S, Q);
        curTime += (double)(clock() - beg) / CLOCKS_PER_SEC;
    } while (curTime < maxTime);

    return {mean, stdDev, ConfidenceInterval(mean, halfDelta), N, curTime};
}

void ImportanceSampling::setSeed(const std::seed_seq &seed) {
    generator.setSeed(seed);
}

void ImportanceSampling::sample(size_t step, size_t& N, double& S, double& Q) {
    const PiecewiseLinearFunction& f = generator.getPWLFunc();

    for (size_t i = 0; i < step; ++i) {
        double X = generator.generate();
        double Y = g(X) / f(X);

        S += Y;
        Q += Y*Y;
    }

    N += step;

    // multiplication a la fin plutot qu'a chaque iteration dans la boucle
    double tmpS = S * f.A;
    double tmpQ = Q*(f.A * f.A);

    mean = tmpS/N;
    double var = tmpQ/N - mean*mean;
    stdDev = sqrt(var/N);
    halfDelta = 1.96 * stdDev;
}


