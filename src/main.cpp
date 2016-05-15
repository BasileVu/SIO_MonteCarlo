#include <cstdlib>
#include <cmath>
#include <functional>
#include <iostream>
#include <iomanip>
#include <list>

#include "montecarlo/UniformSampling.h"
#include "montecarlo/ImportanceSampling.h"
#include "montecarlo/ControlVariableMethod.h"

using namespace std;

/*
void printSampling(const MonteCarloMethod::Sampling& s) {
    cout << " N. de generations : " << s.N << endl;
    cout << " Aire estimee      : " << s.areaEstimator << endl;
    cout << " IC                : " << s.confidenceInterval << endl;
    cout << " Largeur de l'IC   : " << s.confidenceInterval.delta << endl;
    cout << " Temps d'execution : " << s.elapsedTime << "s" << endl;
    cout << endl;
}*/

const string MAX_WIDTH = "Largeur max IC | ";
const string MIN_TIME = "Temps min [s]  | ";
const string HEADER = "N generations | Aire estimee |         IC         | largeur IC | Temps [s]";

void printSampling(double constraint, const MonteCarloMethod::Sampling& s) {
    cout << setw(14) << constraint << " | ";
    cout << setw(13) << s.N << " | ";
    cout << fixed;
    cout << setw(12) << setprecision(5) << s.areaEstimator << " | ";
    cout << setw(18)  << s.confidenceInterval << " | ";
    cout << setw(10)  << setprecision(5) << s.confidenceInterval.delta << " | ";
    cout << setw(9) << setprecision(3) << s.elapsedTime;
    cout.unsetf(ios_base::floatfield);
    cout << endl;
}

int main () {

    // fonction dont on veut estimer l'aire
    MonteCarloMethod::Func g = [](double x) {
        return (25 + x * (x - 6) * (x - 8) * (x - 14) / 25) * exp(sqrt(1 + cos(x*x / 10)));
    };

    // borne inferieure et superieure
    double a = 0, b = 15;

    // nombre de points a utiliser pour la fonction affine par morceaux
    size_t numPointsPWLFunc = 15;

    // nombre d'étapes a faire dans le cas de generations avec largeur d'IC ou temps limite
    size_t step = 100000;

    // graine utilisee pour les generateurs
    seed_seq seed = {24, 512, 42};

    std::list<double> maxWidths;
    for (double i = 1; i > 0.09; i -= 0.1) {
        maxWidths.push_back(i);
    }

    std::list<double> minTimes;
    for (double i = 1; i < 10; ++i) {
        minTimes.push_back(i);
    }


    {
        UniformSampling us(g, a, b);
        us.setSeed(seed);
        cout << "-- Echantillonage uniforme --" << endl;

        cout << MAX_WIDTH << HEADER << endl;
        for (double maxWidth: maxWidths) {
            printSampling(maxWidth, us.sampleWithMaxDelta(maxWidth, step));
        }
        cout << endl;

        cout << MIN_TIME << HEADER << endl;
        for (double minTime: minTimes) {
            printSampling(minTime, us.sampleWithMinTime(minTime, step));
        }
        cout << endl;
    }

    // creation des points de la fonction affine par morceaux
    Points points = Stats::createPoints(numPointsPWLFunc, g, a, b);

    {
        ImportanceSampling is(g, points.xs, points.ys);
        is.setSeed(seed);

        cout << "-- Echantillonage preferentiel --" << endl;

        cout << MAX_WIDTH << HEADER << endl;
        for (double maxWidth: maxWidths) {
            printSampling(maxWidth, is.sampleWithMaxDelta(maxWidth, step));
        }
        cout << endl;

        cout << MIN_TIME << HEADER << endl;
        for (double minTime: minTimes) {
            printSampling(minTime, is.sampleWithMinTime(minTime, step));
        }
        cout << endl;
    }

    {
        ControlVariable cv(g, a, b, points.xs, points.ys);
        cv.setSeed(seed);
        const size_t M = 10000;

        cout << "-- Echantillonage uniforme avec variable de controle --" << endl;
        cout << "Avec M = " << M << " :" << endl;

        cout << MAX_WIDTH << HEADER << endl;
        for (double maxWidth: maxWidths) {
            printSampling(maxWidth, cv.sampleWithMaxDelta(M, maxWidth, step));
        }
        cout << endl;

        cout << MIN_TIME << HEADER << endl;
        for (double minTime: minTimes) {
            printSampling(minTime, cv.sampleWithMinTime(M, minTime, step));
        }
        cout << endl;
    }


    return EXIT_SUCCESS;
}
