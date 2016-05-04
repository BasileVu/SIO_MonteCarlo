#include <ctime>

#include "ControlVariable.h"

ControlVariable::ControlVariable(const Func& g, double a, double b,
                                 const std::vector<double>& xs, const std::vector<double>& ys)
        :
        MonteCarloMethod(g), a(a), b(b),
        distribution(std::uniform_real_distribution<double>(0, 1)),
        h(xs, ys)

{
    // TODO preconditions
    mu = h.A / (b-a);
}

MonteCarloMethod::Sampling ControlVariable::sampleWithSize(size_t M, size_t N) {

    // phase 1: on genere un "petit" echantillon de taille M
    ControlVariable::ResultFirst firstRes = firstStep(M);

    double SV = firstRes.SV, QV = firstRes.QV, c = firstRes.c;
    size_t step = N - M;
    size_t tmpN = M;

    // phase 2 : on poursuit l'echantillonage jusqu'a la taille de N desiree
    // nombre d'echantillons a generer
    ControlVariable::ResultSecond secondRes = secondStep(step, tmpN, c, SV, QV);
    return createSampling(secondRes.meanV, secondRes.halfDelta, N);
}

MonteCarloMethod::Sampling ControlVariable::sampleWithMaxDelta(size_t M, double maxDelta, size_t step) {

    // phase 1: on genere un "petit" echantillon de taille M
    ControlVariable::ResultFirst firstRes = firstStep(M);

    double SV = firstRes.SV, QV = firstRes.QV, c = firstRes.c;
    size_t N = M;

    // phase 2 : on poursuit l'echantillonage jusqu'a la largeur de l'IC desiree
    ResultSecond secondRes;
    do {
        secondRes = secondStep(step , N, c, SV, QV);
    } while (secondRes.halfDelta * 2 > maxDelta);

    return createSampling(secondRes.meanV, secondRes.halfDelta, N);
}

MonteCarloMethod::Sampling ControlVariable::sampleWithMaxTime(size_t M, double maxTime, size_t step) {

    // phase 1: on genere un "petit" echantillon de taille M
    ControlVariable::ResultFirst firstRes = firstStep(M);

    double SV = firstRes.SV, QV = firstRes.QV, c = firstRes.c;
    size_t N = M;
    double curTime = 0;

    // phase 2 : on poursuit l'echantillonage jusqu'a atteindre le temps maximal desire
    ResultSecond secondRes;
    do {
        clock_t beg = clock();
        secondRes = secondStep(step , N, c, SV, QV);
        curTime += (double)(clock() - beg) / CLOCKS_PER_SEC;
    } while (curTime < maxTime);

    return createSampling(secondRes.meanV, secondRes.halfDelta, N);
}

void ControlVariable::setSeed(const std::seed_seq &seed) {
    std::seed_seq copy = seed;
    generator.seed(copy);
}

ControlVariable::ResultFirst ControlVariable::firstStep(size_t M) {

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

    double covYZ = 0;
    for (size_t i = 0; i < M; ++i) {
        covYZ += (yks[i] - meanY) * (zks[i] - mu);
    }

    covYZ /= M;

    double c = -(covYZ / varZ);

    double SV = 0, QV = 0;
    for (size_t i = 0; i < M; ++i) {
        double V = yks[i] + c * (zks[i] - mu);
        SV += V;
        QV += V * V;
    }

    return {SV, QV, c};
}

ControlVariable::ResultSecond ControlVariable::secondStep(size_t step, size_t& N, double c, double& SV, double& QV) {

    for (size_t i = 0; i < step; ++i) {
        double X = distribution(generator) * (b - a) + a;
        double Y = g(X), Z = h(X);
        double V = Y + c * (Z - mu);
        SV += V;
        QV += V * V;
    }

    N += step;

    double meanV = SV / N;
    double varV = (QV / N) - meanV * meanV;
    double halfDelta = 1.96 * (b-a) * sqrt(varV / N);

    return {meanV, halfDelta};
}

MonteCarloMethod::Sampling ControlVariable::createSampling(double meanV, double halfDelta, size_t N) {
    double areaEstimator = (b-a) * meanV;
    return {areaEstimator, ConfidenceInterval(areaEstimator, halfDelta), N};
}




