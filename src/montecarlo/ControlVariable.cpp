#include <vector>
#include <algorithm>

#include "ControlVariable.h"

ControlVariable::ControlVariable(const MonteCarloMethod::Func& g)
        : MonteCarloMethod(g), distribution(std::uniform_real_distribution<double>(0, 1)) {}

MonteCarloMethod::Sampling ControlVariable::sample(size_t M, size_t N, double a, double b,
                                                   const std::vector<double>& xs, const std::vector<double>& ys) {

    // TODO preconditions

    PiecewiseLinearFunction h(xs, ys);
    double mu = h.A / (b-a);

    // phase 1: on genere un "petit" echantillon de taille M
    ControlVariable::ResultFirst firstRes = firstStep(M, a, b, h, mu);

    double SV = firstRes.SV, QV = firstRes.QV, c = firstRes.c;

    // phase 2 : on poursuit l'echantillonage jusqu'a la taille de N desiree
    size_t step = N - M; // nombre d'echantillons a generer
    size_t tmpN = M;
    ControlVariable::ResultSecond secondRes = secondStep(step, tmpN, a, b, h, mu, c, SV, QV);

    double areaEstimator = (b-a) * secondRes.meanV;

    return {areaEstimator, ConfidenceInterval(areaEstimator, secondRes.halfDelta), N};
}

MonteCarloMethod::Sampling ControlVariable::sample(size_t M, double maxDelta, size_t step, double a, double b,
                                                   const std::vector<double>& xs, const std::vector<double>& ys) {

    PiecewiseLinearFunction h(xs, ys);
    double mu = h.A / (b-a);

    // phase 1: on genere un "petit" echantillon de taille M
    ControlVariable::ResultFirst firstRes = firstStep(M, a, b, h, mu);

    double SV = firstRes.SV, QV = firstRes.QV, c = firstRes.c;
    size_t N = M;

    // phase 2 : on poursuit l'echantillonage jusqu'a la largeur de l'IC desiree
    ResultSecond res;
    do {
        res = secondStep(step, N, a, b, h, mu, c, SV, QV);
    } while (res.halfDelta * 2 > maxDelta);

    double areaEstimator = (b-a) * res.meanV;

    return {areaEstimator, ConfidenceInterval(areaEstimator, res.halfDelta), N};
}

void ControlVariable::setSeed(const std::seed_seq &seed) {
    std::seed_seq copy = seed;
    generator.seed(copy);
}

ControlVariable::ResultFirst ControlVariable::firstStep(size_t M, double a, double b, const PiecewiseLinearFunction& h, double mu) {
    // phase 1: on genere un "petit" echantillon de taille M
    std::vector<double> yks, zks;
    yks.reserve(M), zks.reserve(M);

    for (size_t i = 0; i < M; ++i) {
        double X = distribution(generator) * (b-a) + a;
        yks.push_back(g(X));
        zks.push_back(h(X));
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

    double SV = 0, QV = 0;
    for (size_t i = 0; i < M; ++i) {
        double V = yks[i] + c * (zks[i] - mu);
        SV += V;
        QV += V * V;
    }

    return {SV, QV, c};
}

ControlVariable::ResultSecond ControlVariable::secondStep(size_t step, size_t& N, double a, double b,
                                                         const PiecewiseLinearFunction& h, double mu, double c,
                                                         double& SV, double& QV) {
    for (size_t i = 0; i < step; ++i) {
        double X = distribution(generator) * (b - a) + a;
        double Y = g(X), Z = h(X);
        double V = Y + c * (Z - mu);
        SV += V;
        QV += V * V;
    }

    N += step;

    double meanV = SV / N;
    double varV = (QV / N) - meanV;
    double halfDelta = 1.96 * (b-a) * sqrt(varV / N);

    return {meanV, halfDelta};
}


