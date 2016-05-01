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

    seed_seq seed = {24, 512, 42};

    {
        // ~ 601.989
        UniformSampling us(g);
        us.setSeed(seed);

        clock_t start = clock();
        MonteCarloMethod::Sampling s = us.sample(100000, 0, 15);

        cout << s.areaEstimator << ", " << s.confidenceInterval.toString() << endl;
        cout << (double) (clock() - start) / CLOCKS_PER_SEC << "s" << endl;
    }

    Points points = Stats::createPoints(30, g, 0, 15);

    {
        ImportanceSampling is(g);
        is.setSeed(seed);

        clock_t start = clock();
        MonteCarloMethod::Sampling s = is.sample(100000, points.xs, points.ys);

        cout << s.areaEstimator << ", " << s.confidenceInterval.toString() << endl;
        cout << (double) (clock() - start) / CLOCKS_PER_SEC << " s" << endl;
    }

    {

        ControlVariable cv(g);
        cv.setSeed(seed);

        clock_t start = clock();
        MonteCarloMethod::Sampling s = cv.sample(10000, 1000000, 0, 15, points.xs, points.ys);

        cout << s.areaEstimator << ", " << s.confidenceInterval.toString() << endl;
        cout << (double) (clock() - start) / CLOCKS_PER_SEC << " s" << endl;
    }


    return EXIT_SUCCESS;
}
