#include <ctime>

#include "ImportanceSampling.h"

ImportanceSampling::ImportanceSampling(const std::function<double(double)>& g, const std::vector<double>& xs, const std::vector<double>& ys)
        : MonteCarloMethod(g), generator(xs, ys) {}

MonteCarloMethod::Sampling ImportanceSampling::sampleWithSize(size_t N) {

    double S = 0, Q = 0;
    size_t tmpN = 0;

    Result res = sample(N, tmpN, S, Q);

    return {res.mean, ConfidenceInterval(res.mean, res.halfDelta), N};
}

MonteCarloMethod::Sampling ImportanceSampling::sampleWithMaxDelta(double maxDelta, size_t step) {

    double S = 0, Q = 0;
    size_t N = 0;

    // genere des valeurs tant que la largeur de l'intervalle de confiance est plus grande que "maxDelta"
    Result res;
    do {
        res = sample(step, N, S, Q);
    } while (res.halfDelta * 2 > maxDelta);

    return {res.mean, ConfidenceInterval(res.mean, res.halfDelta), N};
}

MonteCarloMethod::Sampling ImportanceSampling::sampleWithMaxTime(double maxTime, size_t step) {

    double S = 0, Q = 0;
    size_t N = 0;
    double curTime = 0;

    // genere des valeurs tant que la largeur de l'intervalle de confiance est plus grande que "maxDelta"
    Result res;
    do {
        clock_t beg = clock();
        res = sample(step, N, S, Q);
        curTime += (double)(clock() - beg) / CLOCKS_PER_SEC;
    } while (curTime < maxTime);

    return {res.mean, ConfidenceInterval(res.mean, res.halfDelta), N};
}

void ImportanceSampling::setSeed(const std::seed_seq &seed) {
    generator.setSeed(seed);
}

ImportanceSampling::Result ImportanceSampling::sample(size_t step, size_t& N, double& S, double& Q) {
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

    double mean = tmpS/N;
    double var = tmpQ/N - mean*mean;

    double halfDelta = 1.96 * sqrt(var/N);

    return {mean, halfDelta};
}


