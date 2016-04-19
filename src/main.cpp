#include <cstdlib>
#include <cmath>
#include <functional>
#include <iostream>
#include <ctime>

#include "UniformSampling.h"

using namespace std;

int main () {

    std::function<double(double)> g = [](double x) {
        return (25 + x * (x - 6) * (x - 8) * (x - 14) / 25) * exp(sqrt(1 + cos(x*x / 10)));
    };

    UniformSampling us(g);
    clock_t start = clock();
    MonteCarloMethod::Sampling s = us.sample();

    cout << s.areaEstimator << ", " << s.confidenceInterval.toString() << endl;
    cout << (double)(clock() - start) / CLOCKS_PER_SEC << endl;



    return EXIT_SUCCESS;
}
