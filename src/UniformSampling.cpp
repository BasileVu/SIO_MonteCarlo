#include <vector>
#include "UniformSampling.h"


UniformSampling::UniformSampling(const std::function<double(double)>& func)
        : MonteCarloMethod(func) {}

MonteCarloMethod::Sampling UniformSampling::sample() {

    const double a = 0, b = 15;

    const size_t nPoints = 100000000;
    double sum = 0, sumSquares = 0;

    for (size_t i = 0; i < nPoints; ++i) {
        double X = distribution(generator) * (b-a);
        double Y = func(X);

        sum += Y;
        sumSquares += Y*Y;
    }

    double mean = sum/nPoints;
    double var = sumSquares/nPoints - mean*mean;
    double areaEstimator = (b-a) * mean;

    double halfDelta = 1.96 * (b-a) * sqrt(var/nPoints);

    // on retourne l'estimateur de l'aire ainsi que l'intervalle de confiance associee.
    return {areaEstimator, {areaEstimator - halfDelta, areaEstimator + halfDelta, halfDelta * 2}};
}