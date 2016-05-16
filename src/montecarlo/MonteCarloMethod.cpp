#include "MonteCarloMethod.h"

MonteCarloMethod::MonteCarloMethod(const std::function<double(double)>& func) : g(func) {}

void MonteCarloMethod::init() {
    sum = 0;
    sumSquares = 0;
    numGen = 0;

    start = clock();
}
