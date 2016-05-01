#include <vector>
#include <algorithm>
#include "ControlVariable.h"

ControlVariable::ControlVariable(const MonteCarloMethod::Func& g)
        : MonteCarloMethod(g), distribution(std::uniform_real_distribution<double>(0, 1)) {}

MonteCarloMethod::Sampling ControlVariable::sample(size_t M, size_t N, double a, double b, const Func& h, double mu) {

    // TODO preconditions

    // phase 1: on genere un "petit" echantillon de taille M
    std::vector<double> ys, zs;
    ys.reserve(M), zs.reserve(M);

    for (size_t i = 0; i < M; ++i) {
        double X = distribution(generator) * (b-a) + a;
        ys.push_back(g(X));
        zs.push_back(h(X));
    }

    double varZ = 0, meanY = 0;
    for (size_t i = 0; i < M; ++i) {
        double tmp = zs[i] - mu;
        varZ += tmp * tmp;
        meanY += ys[i];
    }

    varZ /= M;
    meanY /= M;

    double varYZ = 0;
    for (size_t i = 0; i < M; ++i) {
        varYZ += (ys[i] - meanY) * (zs[i] - mu);
    }

    varYZ /= M;

    double c = -(varYZ / varZ);

    double sumV = 0, squaresV = 0;
    for (size_t i = 0; i < M; ++i) {
        double V = ys[i] + c * (zs[i] - mu);
        sumV += V;
        squaresV += V * V;
    }

    // phase 2 : on poursuit l'échantillonage jusqu'à la taille de N desiree
    for (size_t i = M; i < N; ++i) {
        double X = distribution(generator) * (b-a) + a;
        double Y = g(X), Z = h(X);
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
