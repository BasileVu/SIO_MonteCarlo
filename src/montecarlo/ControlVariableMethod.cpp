#include <ctime>

#include "ControlVariableMethod.h"

ControlVariable::ControlVariable(const Func& g, double a, double b,
                                 const std::vector<double>& xs, const std::vector<double>& ys)
        :
        MonteCarloMethod(g), a(a), b(b),
        uniformDistr(std::uniform_real_distribution<double>(0, 1)),
        h(xs, ys)

{
    if (a >= b) {
        throw std::invalid_argument("Borne inferieure plus grande ou egale a la borne superieure.");
    }
    mu = h.A / (b-a);
}

void ControlVariable::setSamplingSize(size_t M) {
    this->M = M;
}

MonteCarloMethod::Sampling ControlVariable::sampleWithSize(size_t N) {

    // phase 1 : calcul de la constante 'c'
    computeConstant();

    size_t step = N - M;
    numGen = M;

    // phase 2 : on poursuit l'echantillonage jusqu'a la taille de N desiree
    sample(step);
    return createSampling((double)(clock() - start) / CLOCKS_PER_SEC);
}

MonteCarloMethod::Sampling ControlVariable::sampleWithMaxWidth(double maxWidth, size_t step) {

    // phase 1 : calcul de la constante 'c'
    computeConstant();

    numGen = M;

    // phase 2 : on poursuit l'echantillonage jusqu'a la largeur de l'IC desiree
    do {
        sample(step);
    } while (halfWidth * 2 > maxWidth);

    return createSampling((double)(clock() - start) / CLOCKS_PER_SEC);
}

MonteCarloMethod::Sampling ControlVariable::sampleWithMinTime(double minTime, size_t step) {

    // phase 1 : calcul de la constante 'c'
    computeConstant();

    numGen = M;
    double curTime = 0;

    // phase 2 : on poursuit l'echantillonage jusqu'a atteindre le temps minimal desire
    do {
        clock_t beg = clock();
        sample(step);
        curTime += (double)(clock() - beg) / CLOCKS_PER_SEC;
    } while (curTime < minTime);

    return createSampling(curTime);
}

void ControlVariable::setSeed(const std::seed_seq &seed) {
    std::seed_seq copy = seed;
    mtGenerator.seed(copy);
}

void ControlVariable::computeConstant() {

    init();

    std::vector<double> yks, zks;
    yks.reserve(M), zks.reserve(M);

    for (size_t i = 0; i < M; ++i) {
        double X = uniformDistr(mtGenerator) * (b-a) + a; // X ~ U(a,b)
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

void ControlVariable::sample(size_t step) {

    for (size_t i = 0; i < step; ++i) {
        double X = uniformDistr(mtGenerator) * (b - a) + a; // X ~ U(a,b)
        double Y = g(X), Z = h(X);
        double V = Y + c * (Z - mu);
        sum += V;
        sumSquares += V * V;
    }

    numGen += step;

    mean = sum / numGen;
    double var = (sumSquares / numGen) - mean * mean;
    stdDev = (b-a) * sqrt(var / numGen);
    halfWidth = 1.96 * stdDev;
}

MonteCarloMethod::Sampling ControlVariable::createSampling(double timeElapsed) const {
    double areaEstimator = (b-a) * mean;
    return {areaEstimator, stdDev, ConfidenceInterval(areaEstimator, halfWidth), numGen, timeElapsed};
}





