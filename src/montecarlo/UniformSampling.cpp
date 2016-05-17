#include <ctime>

#include <vector>
#include "UniformSampling.h"

UniformSampling::UniformSampling(const MonteCarloMethod::Func& g, double a, double b)
        :
        MonteCarloMethod(g), a(a), b(b),
        uniformDistr(std::uniform_real_distribution<double>(0, 1))
{
    if (b <= a) {
        throw std::invalid_argument("b doit etre plus grand que a");
    }
}

MonteCarloMethod::Sampling UniformSampling::sampleWithSize(uint64_t N) {
    init();

    sample(N);
    return createSampling((double)(clock() - start) / CLOCKS_PER_SEC);
}

MonteCarloMethod::Sampling UniformSampling::sampleWithMaxWidth(double maxWidth, uint64_t step) {
    init();

    do {
        sample(step);
    } while (halfWidth * 2 > maxWidth);

    return createSampling((double)(clock() - start) / CLOCKS_PER_SEC);
}

MonteCarloMethod::Sampling UniformSampling::sampleWithMinTime(double maxTime, uint64_t step) {
    init();
    double curTime = 0;

    do {
        clock_t beg = clock();
        sample(step);
        curTime += (double)(clock() - beg) / CLOCKS_PER_SEC;
    } while (curTime < maxTime);

    return createSampling(curTime);
}

void UniformSampling::setSeed(const std::seed_seq &seed) {
    std::seed_seq copy = seed;
    mtGenerator.seed(copy);
}

void UniformSampling::sample(uint64_t step) {
    for (uint64_t i = 0; i < step; ++i) {
        double X = uniformDistr(mtGenerator) * (b - a) + a; // X ~ U(a,b)
        double Y = g(X);

        sum += Y;
        sumSquares += Y * Y;
    }

    numGen += step;

    mean = sum / numGen;
    double var = (sumSquares / numGen) - mean * mean;
    stdDev = (b-a) * sqrt(var / numGen);
    halfWidth = 1.96 * (b - a) * sqrt(var / numGen);
}

MonteCarloMethod::Sampling UniformSampling::createSampling(double timeElapsed) const {
    double areaEstimator = (b-a) * mean;
    return {areaEstimator, stdDev, ConfidenceInterval(areaEstimator, halfWidth), numGen, timeElapsed};
}