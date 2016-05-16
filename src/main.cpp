#include <cstdlib>
#include <cmath>
#include <functional>
#include <iostream>
#include <iomanip>
#include <list>
#include <fstream>

#include "montecarlo/UniformSampling.h"
#include "montecarlo/ImportanceSampling.h"
#include "montecarlo/ControlVariableMethod.h"

using namespace std;

#define EXPORT_CSV true
const string CSV_FILE = "result.csv";
const char CSV_SEPARATOR = ';';
const string CSV_HEADER = "N generations;Aire estimee;IC;IC inf;IC sup;largeur IC;Temps [s]";

const string MAX_WIDTH = "Largeur max IC";
const string MIN_TIME = "Temps min [s]";
const string HEADER = "N generations | Aire estimee |         IC         | largeur IC | Temps [s]";


void exportSampling(ofstream& ofs, const MonteCarloMethod::Sampling& s) {
    ofs << s.N << CSV_SEPARATOR;
    ofs << s.areaEstimator << CSV_SEPARATOR;
    ofs << s.confidenceInterval << CSV_SEPARATOR;
    ofs << s.confidenceInterval.lower << CSV_SEPARATOR;
    ofs << s.confidenceInterval.upper << CSV_SEPARATOR;
    ofs << s.confidenceInterval.delta << CSV_SEPARATOR;
    ofs << s.elapsedTime << endl;
}

void printSampling(const MonteCarloMethod::Sampling& s) {
    cout << setw(13) << s.N << " | ";
    cout << fixed;
    cout << setw(12) << setprecision(5) << s.areaEstimator << " | ";
    cout << setw(18)  << s.confidenceInterval << " | ";
    cout << setw(10)  << setprecision(5) << s.confidenceInterval.delta << " | ";
    cout << setw(9) << setprecision(3) << s.elapsedTime;
    cout.unsetf(ios_base::floatfield);
    cout << endl;

    if (EXPORT_CSV) {
        ofstream ofs(CSV_FILE, ios_base::app);
        exportSampling(ofs, s);
        ofs.close();
    }
}

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

    if (EXPORT_CSV) {
        ofstream ofs(CSV_FILE, ios_base::app);
        ofs << constraint << CSV_SEPARATOR;
        exportSampling(ofs, s);
        ofs.close();
    }
}

int main () {

    // fonction dont on veut estimer l'aire
    MonteCarloMethod::Func g = [](double x) {
        return (25 + x * (x - 6) * (x - 8) * (x - 14) / 25) * exp(sqrt(1 + cos(x*x / 10)));
    };

    // borne inferieure et superieure
    double a = 0, b = 15;

    // creation des points de la fonction affine par morceaux
    Points points = Stats::createPoints(15, g, a, b);

    // nombre d'étapes a faire dans le cas de generations avec largeur d'IC ou temps limite
    size_t step = 100000;

    // graine utilisee pour les generateurs
    seed_seq seed = {24, 512, 42};

    list<double> maxWidths;
    for (double i = 1; i > 0.09; i -= 0.1) {
        maxWidths.push_back(i);
    }
    for (double i = 0.09; i > 0.04; i -= 0.01) {
        maxWidths.push_back(i);
    }

    list<double> minTimes;
    for (double i = 1; i < 30; ++i) {
        minTimes.push_back(i);
    }

    // cree un nouveau fichier / ecrase l'ancien s'il existe
    if (EXPORT_CSV) {
        ofstream ofs(CSV_FILE);
        ofs.close();
    }


    cout << "-----------------------------------------------------" << endl;
    cout << "| Test de l'implementation des differentes methodes |" << endl;
    cout << "-----------------------------------------------------" << endl << endl;
    {
        UniformSampling us(g, a, b);
        ImportanceSampling is(g, points.xs, points.ys);
        ControlVariable cv(g, a, b, points.xs, points.ys);

        us.setSeed(seed);
        is.setSeed(seed);
        cv.setSeed(seed);

        const size_t N = 100000, M = 10000;

        cout << "-- Echantillonage uniforme --" << endl;
        cout << HEADER << endl;
        if (EXPORT_CSV) {
            ofstream ofs(CSV_FILE, ios_base::app);
            ofs << CSV_HEADER << endl;
            ofs.close();
        }

        for (size_t i = N; i < N * 1000; i *= 10) {
            printSampling(us.sampleWithSize(i));
        }
        cout << endl;

        cout << "-- Echantillonage preferentiel --" << endl;
        cout << HEADER << endl;
        if (EXPORT_CSV) {
            ofstream ofs(CSV_FILE, ios_base::app);
            ofs << CSV_HEADER << endl;
            ofs.close();
        }

        for (size_t i = N; i < N * 1000; i *= 10) {
            printSampling(is.sampleWithSize(i));
        }
        cout << endl;

        cout << "-- Echantillonage uniforme avec variable de controle --" << endl;
        cout << HEADER << endl;
        if (EXPORT_CSV) {
            ofstream ofs(CSV_FILE, ios_base::app);
            ofs << CSV_HEADER << endl;
            ofs.close();
        }
        for (size_t i = N; i < N * 1000; i *= 10) {
            printSampling(cv.sampleWithSize(M, i));
        }
        cout << endl << endl;
    }

    cout << "---------------------------------------------------------------------" << endl;
    cout << "| Echantillonages en utillisant differentes methodes et contraintes |" << endl;
    cout << "---------------------------------------------------------------------" << endl << endl;

    {
        UniformSampling us(g, a, b);
        us.setSeed(seed);
        cout << "-- Echantillonage uniforme --" << endl;

        cout << MAX_WIDTH << " | " << HEADER << endl;
        if (EXPORT_CSV) {
            ofstream ofs(CSV_FILE, ios_base::app);
            ofs << MAX_WIDTH << CSV_SEPARATOR << CSV_HEADER << endl;
            ofs.close();
        }
        for (double maxWidth: maxWidths) {
            printSampling(maxWidth, us.sampleWithMaxDelta(maxWidth, step));
        }
        cout << endl;

        cout << MIN_TIME << "  | " << HEADER << endl;
        if (EXPORT_CSV) {
            ofstream ofs(CSV_FILE, ios_base::app);
            ofs << MIN_TIME << CSV_SEPARATOR << CSV_HEADER << endl;
            ofs.close();
        }
        for (double minTime: minTimes) {
            printSampling(minTime, us.sampleWithMinTime(minTime, step));
        }
        cout << endl;
    }

    {
        ImportanceSampling is(g, points.xs, points.ys);
        is.setSeed(seed);

        cout << "-- Echantillonage preferentiel --" << endl;

        cout << MAX_WIDTH << " | " << HEADER << endl;
        if (EXPORT_CSV) {
            ofstream ofs(CSV_FILE, ios_base::app);
            ofs << MAX_WIDTH << CSV_SEPARATOR << CSV_HEADER << endl;
            ofs.close();
        }
        for (double maxWidth: maxWidths) {
            printSampling(maxWidth, is.sampleWithMaxDelta(maxWidth, step));
        }
        cout << endl;

        cout << MIN_TIME << "  | " << HEADER << endl;
        if (EXPORT_CSV) {
            ofstream ofs(CSV_FILE, ios_base::app);
            ofs << MIN_TIME << CSV_SEPARATOR << CSV_HEADER << endl;
            ofs.close();
        }
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

        cout << MAX_WIDTH << " | " << HEADER << endl;
        if (EXPORT_CSV) {
            ofstream ofs(CSV_FILE, ios_base::app);
            ofs << MAX_WIDTH << CSV_SEPARATOR << CSV_HEADER << endl;
            ofs.close();
        }
        for (double maxWidth: maxWidths) {
            printSampling(maxWidth, cv.sampleWithMaxDelta(M, maxWidth, step));
        }
        cout << endl;

        cout << MIN_TIME << "  | " << HEADER << endl;
        if (EXPORT_CSV) {
            ofstream ofs(CSV_FILE, ios_base::app);
            ofs << MIN_TIME << CSV_SEPARATOR << CSV_HEADER << endl;
            ofs.close();
        }
        for (double minTime: minTimes) {
            printSampling(minTime, cv.sampleWithMinTime(M, minTime, step));
        }
        cout << endl;
    }


    return EXIT_SUCCESS;
}
