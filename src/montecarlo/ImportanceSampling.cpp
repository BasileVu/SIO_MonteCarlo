#include "ImportanceSampling.h"
#include "../generators/RandomValueGenerator.h"

ImportanceSampling::ImportanceSampling(const std::function<double(double)>& g)
        : MonteCarloMethod(g) {}

MonteCarloMethod::Sampling ImportanceSampling::sample(size_t N, const std::vector<double> xs, const std::vector<double> ys) {

    double S = 0, Q = 0;
    size_t tmpN = N;

    InverseFunctions inv(xs, ys);
    inv.setSeed(seed);
    Result res = sample(N, tmpN, inv, S, Q);

    return {res.mean, ConfidenceInterval(res.mean, res.halfDelta), N};
}

MonteCarloMethod::Sampling ImportanceSampling::sample(double maxDelta, size_t step,
                                    const std::vector<double> xs, const std::vector<double> ys) {

    double S = 0, Q = 0;
    size_t N = 0;

    InverseFunctions inv(xs, ys);
    inv.setSeed(seed);

    // genere des valeurs tant que la largeur de l'intervalle de confiance est plus grande que "maxDelta"
    Result res;
    do {
        res = sample(step, N, inv, S, Q);
    } while (res.halfDelta * 2 > maxDelta);

    return {res.mean, ConfidenceInterval(res.mean, res.halfDelta), N};
}

void ImportanceSampling::setSeed(const std::seed_seq &seed) {
    this->seed = seed;
}

ImportanceSampling::Result ImportanceSampling::sample(size_t step, size_t& N, RandomValueGenerator& rvg, double& S, double& Q) {
    const PiecewiseLinearFunction& f = rvg.getPWLFunc();

    for (size_t i = 0; i < step; ++i) {
        double X = rvg.generate();
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


