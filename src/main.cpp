#include <cstdlib>
#include <cmath>
#include <functional>
#include <iostream>
#include <ctime>

#include "montecarlo/UniformSampling.h"
#include "montecarlo/ImportanceSampling.h"
#include "montecarlo/ControlVariable.h"

using namespace std;

int main () {

    MonteCarloMethod::Func g = [](double x) {
        //return sqrt(x + cos(x*x));
        return (25 + x * (x - 6) * (x - 8) * (x - 14) / 25) * exp(sqrt(1 + cos(x*x / 10)));
    };

    double a = 0, b = 15;
    size_t M = 10000, N = 1000000;
    size_t numPointsPWL = 30;

    seed_seq seed = {24, 512, 42};

    {
        UniformSampling us(g);
        us.setSeed(seed);

        clock_t start = clock();
        MonteCarloMethod::Sampling s = us.sample(N, a, b);

        cout << s.areaEstimator << ", " << s.confidenceInterval.toString() << endl;
        cout << (double) (clock() - start) / CLOCKS_PER_SEC << "s" << endl;
    }

    Points points = Stats::createPoints(numPointsPWL, g, a, b);

    {
        ImportanceSampling is(g);
        is.setSeed(seed);

        clock_t start = clock();
        MonteCarloMethod::Sampling s = is.sample(N, points.xs, points.ys);

        cout << s.areaEstimator << ", " << s.confidenceInterval.toString() << endl;
        cout << (double) (clock() - start) / CLOCKS_PER_SEC << " s" << endl;
    }

    {
        ControlVariable cv(g);
        cv.setSeed(seed);

        clock_t start = clock();
        MonteCarloMethod::Sampling s = cv.sample(M, N, a, b, points.xs, points.ys);

        cout << s.areaEstimator << ", " << s.confidenceInterval.toString() << endl;
        cout << (double) (clock() - start) / CLOCKS_PER_SEC << " s" << endl;
    }


    return EXIT_SUCCESS;
}
