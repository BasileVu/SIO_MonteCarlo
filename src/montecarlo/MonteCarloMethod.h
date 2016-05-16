#ifndef MONTECARLOMETHOD_H
#define MONTECARLOMETHOD_H

#include <random>
#include <functional>
#include <ctime>

#include "../utility/Stats.h"

/**
 * Represente une methode de Monte-Carlo (dans notre cas, utilisee afin de calculer une integrale en estimant son aire).
 */
class MonteCarloMethod {
public:
    // une fonction prenant un double et retournant un double
    typedef std::function<double(double)> Func;

protected:
    const Func& g;      // la fonciton dont on veut evaluer l'aire

    double mean;        // la moyenne des valeurs generees
    double stdDev;      // l'estimateur de l'ecart-type de l'estimateur de l'aire
    double halfWidth;   // la demi-largeur de l'IC pour l'estimateur de l'aire

    double sum;         // la somme des valeurs
    double sumSquares;  // la somme des carres des valeurs
    size_t numGen;      // la taille de l'echantillon (nombre de valeurs generees)

    clock_t start;      // utile pour la mesure du temps requis pour generer un echantillon

public:
    /**
     * Represente le resultat d'un echantillonnage.
     */
    struct Sampling {
        double areaEstimator;                   // aire estimee
        double stdDevEstimator;                 // estimateur de l'ecart-type de l'aire estimee
        ConfidenceInterval confidenceInterval;  // intervalle de confiance a 95%
        size_t N;                               // taille de l'echantillon
        double elapsedTime;                     // temps pour creer la totalite de l'echantillon
    };

    /*
     * Prepare la methode.
     *
     * @param g La fonction dont on veut estimer l'aire.
     */
    MonteCarloMethod(const Func& g);

    /**
     * Initialise la graine du generateur utilise pour la methode.
     */
    virtual void setSeed(const std::seed_seq& seed) = 0;

protected:
    /**
     * Initialise les differents champs. Doit etre appelee au debut de chaque etape d'echantillonage.
     */
    void init();
};

#endif // MONTECARLOMETHOD_H
