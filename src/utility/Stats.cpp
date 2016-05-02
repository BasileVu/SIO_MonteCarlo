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

double Stats::expectedValue(const PiecewiseLinearFunction& f) {
    double res = 0;
    for (const Piece& p : f.pieces) {
        double x0 = p.x0, x1 = p.x1;
        double y0 = p.y0, y1 = p.y1;

        if (y0 + y1 > 0) {
            double eVal = (y0 * (2 * x0 + x1) + y1 * (x0 + 2 * x1)) / (3 * (y0 + y1));
            res += eVal * p.A_k;
        }
    }
    return res / f.A;
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

    if (numPoints < 2) {
        throw std::invalid_argument("Le nombre de points doit etre au moins egal a 2.");
    }

    if (a > b) {
        throw std::invalid_argument("Borne minimale plus grande que borne maximale.");
    }

    std::vector<double> xs, ys;
    xs.reserve(numPoints); ys.reserve(numPoints);

    double pieceWidth = (b - a) / (numPoints - 1);

    for (size_t i = 0; i < numPoints; ++i) {
        double x = a + pieceWidth * i;
        xs.push_back(x);
        ys.push_back(func(x));
    }

    return {xs, ys};
}