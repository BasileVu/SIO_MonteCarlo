#include <sstream>

#include "Stats.h"

std::string ConfidenceInterval::toString() const {
    std::stringstream ss;
    ss << "[" << lower << "," << upper << "]";
    return ss.str();
}

double Stats::sampleVar(const std::vector<double>& values) {
    double sum = 0;
    double m = mean(values);

    for (const double& d : values) {
        double item = d - m;
        sum += item * item;
    }
    return sum / (values.size() - 1);
}

double Stats::sampleStdDev(const std::vector<double>& values) {
    return std::sqrt(sampleVar(values));
}


double Stats::mean(const std::vector<double>& values) {
    return std::accumulate(values.begin(), values.end(), 0.0) / values.size();
}

ConfidenceInterval Stats::confidenceInterval(const std::vector<double>& values, double quantile) {
    double m = mean(values);
    double halfDelta = quantile * (sampleStdDev(values) / std::sqrt(values.size()));
    return ConfidenceInterval {m - halfDelta, m + halfDelta, halfDelta*2};
};

Points Stats::createPoints(size_t numPoints, const std::function<double(double)>& func, double a, double b) {

    if (numPoints == 0) {
        throw std::invalid_argument("Le nombre de points doit etre plus grand que 0.");
    }

    if (a > b) {
        throw std::invalid_argument("a est plus grand que b !");
    }

    std::vector<double> xs, ys;
    xs.reserve(numPoints); ys.reserve(numPoints);

    double pieceWidth = (b - a)/numPoints;

    for (size_t i = 0; i < numPoints; ++i) {
        double x = a + pieceWidth * i;
        xs.push_back(x);
        ys.push_back(func(x));
    }

    return {xs, ys};
}