#include <vector>
#include <algorithm>

#include "ControlVariable.h"

ControlVariable::ControlVariable(const MonteCarloMethod::Func& g)
        : MonteCarloMethod(g), distribution(std::uniform_real_distribution<double>(0, 1)) {}

MonteCarloMethod::Sampling ControlVariable::sample(size_t M, size_t N, double a, double b,
                                                   const std::vector<double>& xs, const std::vector<double>& ys) {

    // TODO preconditions

    PiecewiseLinearFunction h(xs, ys);
    double mu = Stats::expectedValue(h);

    // phase 1: on genere un "petit" echantillon de taille M
    std::vector<double> yks, zks;
    yks.reserve(M), zks.reserve(M);

    for (size_t i = 0; i < M; ++i) {
        double X = distribution(generator) * (b-a) + a;
        yks.push_back(g(X));
        zks.push_back(h(X) / h.A);
    }

    double varZ = 0, meanY = 0;
    for (size_t i = 0; i < M; ++i) {
        double tmp = zks[i] - mu;
        varZ += tmp * tmp;
        meanY += yks[i];
    }

    varZ /= M;
    meanY /= M;

    double varYZ = 0;
    for (size_t i = 0; i < M; ++i) {
        varYZ += (yks[i] - meanY) * (zks[i] - mu);
    }

    varYZ /= M;

    double c = -(varYZ / varZ);

    double sumV = 0, squaresV = 0;
    for (size_t i = 0; i < M; ++i) {
        double V = yks[i] + c * (zks[i] - mu);
        sumV += V;
        squaresV += V * V;
    }

    // phase 2 : on poursuit l'echantillonage jusqu'a la taille de N desiree
    for (size_t i = M; i < N; ++i) {
        double X = distribution(generator) * (b-a) + a;
        double Y = g(X), Z = h(X) / h.A;
        double V = Y + c * (Z - mu);
        sumV += V;
        squaresV += V * V;
    }

    double meanV = sumV / N;
    double varV = (squaresV / N) - meanV;
    double areaEstimator = (b-a) * meanV;
    double halfDelta = 1.96 * (b-a) * sqrt(varV / N);

    return {areaEstimator, {areaEstimator - halfDelta, areaEstimator + halfDelta, halfDelta * 2}};
}

void ControlVariable::setSeed(const std::seed_seq &seed) {
    std::seed_seq copy = seed;
    generator.seed(copy);
}
