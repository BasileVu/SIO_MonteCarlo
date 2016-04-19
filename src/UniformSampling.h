#ifndef UNIFORMSAMPLING_H
#define UNIFORMSAMPLING_H

#include "MonteCarloMethod.h"

class UniformSampling : public MonteCarloMethod {
public:
    UniformSampling(const std::function<double(double)>& func);

    Sampling sample();
};

#endif // UNIFORMSAMPLING_H
