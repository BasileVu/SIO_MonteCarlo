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
    firstStep(M);

    size_t step = N - M;
    numGen = M;

    // phase 2 : on poursuit l'echantillonage jusqu'a la taille de N desiree
    // nombre d'echantillons a generer
    secondStep(step);
    return createSampling((double)(clock() - start) / CLOCKS_PER_SEC);
}

MonteCarloMethod::Sampling ControlVariable::sampleWithMaxDelta(size_t M, double maxDelta, size_t step) {

    // phase 1: on genere un "petit" echantillon de taille M
    firstStep(M);

    numGen = M;

    // phase 2 : on poursuit l'echantillonage jusqu'a la largeur de l'IC desiree
    do {
        secondStep(step);
    } while (halfDelta * 2 > maxDelta);

    return createSampling((double)(clock() - start) / CLOCKS_PER_SEC);
}

MonteCarloMethod::Sampling ControlVariable::sampleWithMinTime(size_t M, double maxTime, size_t step) {

    // phase 1: on genere un "petit" echantillon de taille M
    firstStep(M);

    numGen = M;
    double curTime = 0;

    // phase 2 : on poursuit l'echantillonage jusqu'a atteindre le temps maximal desire
    do {
        clock_t beg = clock();
        secondStep(step);
        curTime += (double)(clock() - beg) / CLOCKS_PER_SEC;
    } while (curTime < maxTime);

    return createSampling(curTime);
}

void ControlVariable::setSeed(const std::seed_seq &seed) {
    std::seed_seq copy = seed;
    generator.seed(copy);
}

void ControlVariable::firstStep(size_t M) {

    init();

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

    c = -(covYZ / varZ);

    for (size_t i = 0; i < M; ++i) {
        double V = yks[i] + c * (zks[i] - mu);
        sum += V;
        sumSquares += V * V;
    }
}

void ControlVariable::secondStep(size_t step) {

    for (size_t i = 0; i < step; ++i) {
        double X = distribution(generator) * (b - a) + a;
        double Y = g(X), Z = h(X);
        double V = Y + c * (Z - mu);
        sum += V;
        sumSquares += V * V;
    }

    numGen += step;

    mean = sum / numGen;
    double varV = (sumSquares / numGen) - mean * mean;
    stdDev = (b-a) * sqrt(varV / numGen);
    halfDelta = 1.96 * stdDev;
}

MonteCarloMethod::Sampling ControlVariable::createSampling(double timeElapsed) const {
    double areaEstimator = (b-a) * mean;
    return {areaEstimator, stdDev, ConfidenceInterval(areaEstimator, halfDelta), numGen, timeElapsed};
}




