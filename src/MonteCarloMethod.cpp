#include "MonteCarloMethod.h"

MonteCarloMethod::MonteCarloMethod(const std::function<double(double)> &func)
        : func(func), distribution(std::uniform_real_distribution<double>(0, 1)) {}
