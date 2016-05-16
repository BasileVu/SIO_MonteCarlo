#include <ctime>

#include "ControlVariableMethod.h"

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

    clock_t start = clock();

    // phase 1: on genere un "petit" echantillon de taille M
    ControlVariable::ResultFirst firstRes = firstStep(M);

    double SV = firstRes.SV, QV = firstRes.QV, c = firstRes.c;
    size_t step = N - M;
    size_t tmpN = M;

    // phase 2 : on poursuit l'echantillonage jusqu'a la taille de N desiree
    // nombre d'echantillons a generer
    secondStep(step, tmpN, c, SV, QV);
    return createSampling(N, (double)(clock() - start) / CLOCKS_PER_SEC);
}

MonteCarloMethod::Sampling ControlVariable::sampleWithMaxDelta(size_t M, double maxDelta, size_t step) {

    clock_t start = clock();

    // phase 1: on genere un "petit" echantillon de taille M
    ControlVariable::ResultFirst firstRes = firstStep(M);

    double SV = firstRes.SV, QV = firstRes.QV, c = firstRes.c;
    size_t N = M;

    // phase 2 : on poursuit l'echantillonage jusqu'a la largeur de l'IC desiree
    do {
        secondStep(step , N, c, SV, QV);
    } while (halfDelta * 2 > maxDelta);

    return createSampling(N, (double)(clock() - start) / CLOCKS_PER_SEC);
}

MonteCarloMethod::Sampling ControlVariable::sampleWithMinTime(size_t M, double maxTime, size_t step) {

    // phase 1: on genere un "petit" echantillon de taille M
    ControlVariable::ResultFirst firstRes = firstStep(M);

    double SV = firstRes.SV, QV = firstRes.QV, c = firstRes.c;
    size_t N = M;
    double curTime = 0;

    // phase 2 : on poursuit l'echantillonage jusqu'a atteindre le temps maximal desire
    do {
        clock_t beg = clock();
        secondStep(step , N, c, SV, QV);
        curTime += (double)(clock() - beg) / CLOCKS_PER_SEC;
    } while (curTime < maxTime);

    return createSampling(N, curTime);
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

void ControlVariable::secondStep(size_t step, size_t& N, double c, double& SV, double& QV) {

    for (size_t i = 0; i < step; ++i) {
        double X = distribution(generator) * (b - a) + a;
        double Y = g(X), Z = h(X);
        double V = Y + c * (Z - mu);
        SV += V;
        QV += V * V;
    }

    N += step;

    mean = SV / N;
    double varV = (QV / N) - mean * mean;
    stdDev = (b-a) * sqrt(varV / N);
    halfDelta = 1.96 * stdDev;
}

MonteCarloMethod::Sampling ControlVariable::createSampling(size_t N, double timeElapsed) {
    double areaEstimator = (b-a) * mean;
    return {areaEstimator, stdDev, ConfidenceInterval(areaEstimator, halfDelta), N, timeElapsed};
}




