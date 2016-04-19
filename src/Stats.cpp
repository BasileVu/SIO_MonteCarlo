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
