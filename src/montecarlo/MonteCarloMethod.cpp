#include "MonteCarloMethod.h"

MonteCarloMethod::MonteCarloMethod(const std::function<double(double)>& func) : g(func) {}