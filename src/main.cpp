#include <cstdlib>
#include <cmath>
#include <functional>
#include <iostream>
#include <ctime>
#include <string>

#include "montecarlo/UniformSampling.h"
#include "montecarlo/ImportanceSampling.h"
#include "montecarlo/ControlVariableMethod.h"

using namespace std;

void printSampling(const MonteCarloMethod::Sampling& s, const string& name, double timeTaken) {
    cout << "-- " << name << " --" << endl;
    cout << " N. de generations : " << s.N << endl;
    cout << " Aire estimee      : " << s.areaEstimator << endl;
    cout << " IC                : " << s.confidenceInterval << endl;
    cout << " Temps d'execution : " << timeTaken << "s" << endl;
    cout << endl;
}

int main () {

    // fonction dont on veut estimer l'aire
    MonteCarloMethod::Func g = [](double x) {
        return (25 + x * (x - 6) * (x - 8) * (x - 14) / 25) * exp(sqrt(1 + cos(x*x / 10)));
    };

    // borne inferieure et superieure
    double a = 0, b = 15;

    // M : nombre de valeurs dans l'echantillon, N : nombre de valeurs totales a generer
    size_t M = 10000, N = 1000000;

    // nombre de points a utiliser pour la fonction affine par morceaux
    size_t numPointsPWLFunc = 15;

    // largeur max de l'IC à 95%
    double deltaMax = 1;

    // temps maximum alloue pour chaque methode (en secondes)
    double maxTime = 2;

    size_t step = 100000;

    // graine utilisee pour les generateurs
    seed_seq seed = {24, 512, 42};

    cout << "Largeur max de l'IC : " << deltaMax << endl;
    cout << endl;

    {
        UniformSampling us(g, a, b);
        us.setSeed(seed);

        clock_t start = clock();
        MonteCarloMethod::Sampling s = us.sampleWithMaxTime(maxTime, step);
        double timeTaken = (double)(clock() - start) / CLOCKS_PER_SEC;

        printSampling(s, "Echantillonage uniforme", timeTaken);

    }

    // creation des points de la fonction affine par morceaux
    Points points = Stats::createPoints(numPointsPWLFunc, g, a, b);

    {
        ImportanceSampling is(g, points.xs, points.ys);
        is.setSeed(seed);

        clock_t start = clock();
        MonteCarloMethod::Sampling s = is.sampleWithMaxTime(maxTime, step);
        double timeTaken = (double)(clock() - start) / CLOCKS_PER_SEC;

        printSampling(s, "Echantillonage preferentiel", timeTaken);
    }

    {
        ControlVariable cv(g, a, b, points.xs, points.ys);
        cv.setSeed(seed);

        clock_t start = clock();
        MonteCarloMethod::Sampling s = cv.sampleWithMaxTime(M, maxTime, step);
        double timeTaken = (double)(clock() - start) / CLOCKS_PER_SEC;

        printSampling(s, "Methode de la variable de controle", timeTaken);
    }


    return EXIT_SUCCESS;
}
