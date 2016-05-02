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
    ControlVariable::Result res = firstStep(M, a, b, h, mu);

    double SV = res.SV, QV = res.QV, c = res.c;

    // phase 2 : on poursuit l'echantillonage jusqu'a la taille de N desiree
    for (size_t i = M; i < N; ++i) {
        double X = distribution(generator) * (b-a) + a;
        double Y = g(X), Z = h(X) / h.A;
        double V = Y + c * (Z - mu);
        SV += V;
        QV += V * V;
    }

    double meanV = SV / N;
    double varV = (SV / N) - meanV;
    double halfDelta = 1.96 * (b-a) * sqrt(varV / N);

    double areaEstimator = (b-a) * meanV;

    return {areaEstimator, ConfidenceInterval(areaEstimator, halfDelta), N};
}

MonteCarloMethod::Sampling ControlVariable::sample(size_t M, double maxDelta, size_t step, double a, double b,
                                                   const std::vector<double>& xs, const std::vector<double>& ys) {

    PiecewiseLinearFunction h(xs, ys);
    double mu = Stats::expectedValue(h);

    // phase 1: on genere un "petit" echantillon de taille M
    ControlVariable::Result res = firstStep(M, a, b, h, mu);

    double SV = res.SV, QV = res.QV, c = res.c;
    size_t N = M;

    // phase 2 : on poursuit l'echantillonage jusqu'a lalageur de l'IC desiree
    double halfDelta = 0, meanV = 0;
    do {
        for (size_t i = 0; i < step; ++i) {
            double X = distribution(generator) * (b - a) + a;
            double Y = g(X), Z = h(X) / h.A;
            double V = Y + c * (Z - mu);
            SV += V;
            QV += V * V;
        }

        N += step;

        meanV = SV / N;
        double varV = (QV / N) - meanV;
        halfDelta = 1.96 * (b-a) * sqrt(varV / N);

    } while (halfDelta * 2 > maxDelta);

    double areaEstimator = (b-a) * meanV;

    return {areaEstimator, ConfidenceInterval(areaEstimator, halfDelta), N};
}

void ControlVariable::setSeed(const std::seed_seq &seed) {
    std::seed_seq copy = seed;
    generator.seed(copy);
}

ControlVariable::Result ControlVariable::firstStep(size_t M, double a, double b, const PiecewiseLinearFunction& h, double mu) {
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

    double SV = 0, QV = 0;
    for (size_t i = 0; i < M; ++i) {
        double V = yks[i] + c * (zks[i] - mu);
        SV += V;
        QV += V * V;
    }

    return {SV, QV, c};
}


