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
const string CSV_HEADER = "N generations;Aire estimee;IC;RC(N) * ET;largeur IC;Temps [s];IC inf;IC sup";

const string MAX_WIDTH = "Largeur max IC";
const string MIN_TIME = "Temps min [s]";
const string HEADER = "N generations | Aire estimee |      IC a 95%      | RC(N) * ET | largeur IC | Temps [s]";


void exportSampling(ofstream& ofs, const MonteCarloMethod::Sampling& s) {
    ofs << s.N << CSV_SEPARATOR;
    ofs << fixed;
    ofs << setprecision(5) << s.areaEstimator << CSV_SEPARATOR;
    ofs << s.confidenceInterval << CSV_SEPARATOR;
    ofs << (sqrt(s.N) * s.stdDevEstimator) << CSV_SEPARATOR;
    ofs << s.confidenceInterval.width << CSV_SEPARATOR;
    ofs << setprecision(3) << s.elapsedTime << CSV_SEPARATOR;
    ofs << s.confidenceInterval.lower << CSV_SEPARATOR;
    ofs << s.confidenceInterval.upper << CSV_SEPARATOR;
    ofs << endl;
}

void printSampling(const MonteCarloMethod::Sampling& s) {
    cout << setw(13) << s.N << " | ";
    cout << fixed;
    cout << setw(12) << setprecision(5) << s.areaEstimator << " | ";
    cout << setw(18)  << s.confidenceInterval << " | ";
    cout << setw(10)  << (sqrt(s.N) * s.stdDevEstimator) << " | ";
    cout << setw(10)  << s.confidenceInterval.width << " | ";
    cout << setw(9) << setprecision(3) << s.elapsedTime;
    cout.unsetf(ios_base::floatfield);
    cout << endl;
}

void printExportSampling(const MonteCarloMethod::Sampling& s) {
    printSampling(s);

    if (EXPORT_CSV) {
        ofstream ofs(CSV_FILE, ios_base::app);
        exportSampling(ofs, s);
        ofs.close();
    }
}

void printExportSampling(double constraint, const MonteCarloMethod::Sampling& s) {
    cout << setw(14) << constraint << " | ";
    printSampling(s);

    if (EXPORT_CSV) {
        ofstream ofs(CSV_FILE, ios_base::app);
        ofs << constraint << CSV_SEPARATOR;
        exportSampling(ofs, s);
        ofs.close();
    }
}

void runImplementationTest(MonteCarloMethod& m) {
    cout << HEADER << endl;
    if (EXPORT_CSV) {
        ofstream ofs(CSV_FILE, ios_base::app);
        ofs << CSV_HEADER << endl;
        ofs.close();
    }

    for (size_t i = 100000; i <= 10000000; i *= 10) {
        printExportSampling(m.sampleWithSize(i));
    }
    cout << endl;
}

void runTests(MonteCarloMethod& m, const list<double>& maxWidths, const list<double>& minTimes, size_t step) {

    cout << MAX_WIDTH << " | " << HEADER << endl;
    if (EXPORT_CSV) {
        ofstream ofs(CSV_FILE, ios_base::app);
        ofs << MAX_WIDTH << CSV_SEPARATOR << CSV_HEADER << endl;
        ofs.close();
    }
    for (double maxWidth: maxWidths) {
        printExportSampling(maxWidth, m.sampleWithMaxWidth(maxWidth, step));
    }
    cout << endl;

    cout << MIN_TIME << "  | " << HEADER << endl;
    if (EXPORT_CSV) {
        ofstream ofs(CSV_FILE, ios_base::app);
        ofs << MIN_TIME << CSV_SEPARATOR << CSV_HEADER << endl;
        ofs.close();
    }
    for (double minTime: minTimes) {
        printExportSampling(minTime, m.sampleWithMinTime(minTime, step));
    }
    cout << endl;
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

    // nombre d'etapes a faire dans le cas de generations avec une largeur d'IC ou un temps limite
    size_t step = 100000;

    // graine utilisee pour les generateurs
    seed_seq seed = {24, 512, 42};

    // creation des largeurs max
    list<double> maxWidths;
    for (double i = 1; i >= 0.1; i -= 0.1) {
        maxWidths.push_back(i);
    }
    for (double i = 0.09; i >= 0.05; i -= 0.01) {
        maxWidths.push_back(i);
    }

    // creation des temps min
    list<double> minTimes;
    for (double i = 1; i < 30; ++i) {
        minTimes.push_back(i);
    }

    if (EXPORT_CSV) {
        // cree un nouveau fichier / ecrase l'ancien s'il existe
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

        const size_t M = 10000;
        cv.setSamplingSize(M);

        cout << "-- Echantillonage uniforme --" << endl;
        runImplementationTest(us);

        cout << "-- Echantillonage preferentiel --" << endl;
        runImplementationTest(is);

        cout << "-- Echantillonage uniforme avec variable de controle --" << endl;
        runImplementationTest(cv);
    }

    cout << "---------------------------------------------------------------------" << endl;
    cout << "| Echantillonages en utillisant differentes methodes et contraintes |" << endl;
    cout << "---------------------------------------------------------------------" << endl << endl;

    {
        UniformSampling us(g, a, b);
        us.setSeed(seed);
        cout << "-- Echantillonage uniforme --" << endl;
        runTests(us, maxWidths, minTimes, step);
    }

    {
        ImportanceSampling is(g, points.xs, points.ys);
        is.setSeed(seed);

        cout << "-- Echantillonage preferentiel --" << endl;
        runTests(is, maxWidths, minTimes, step);
    }

    {
        ControlVariable cv(g, a, b, points.xs, points.ys);
        cv.setSeed(seed);

        size_t M = 10000;
        cv.setSamplingSize(M);

        cout << "-- Echantillonage uniforme avec variable de controle --" << endl;
        cout << "Avec M = " << M << " :" << endl;
        runTests(cv, maxWidths, minTimes, step);
    }


    return EXIT_SUCCESS;
}
