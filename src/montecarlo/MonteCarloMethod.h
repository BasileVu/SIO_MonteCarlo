#ifndef MONTECARLOMETHOD_H
#define MONTECARLOMETHOD_H

#include <random>
#include <functional>
#include <ctime>
#include <cstdint>

#include "../utility/Stats.h"

/**
 * Represente une methode de Monte-Carlo (dans notre cas, utilisee afin de calculer une integrale en estimant son aire).
 */
class MonteCarloMethod {
public:
    // une fonction prenant un double et retournant un double
    typedef std::function<double(double)> Func;

protected:
    const Func& g;      // la fonciton dont on veut estimer l'aire

    double mean;        // la moyenne des valeurs generees
    double stdDev;      // l'estimateur de l'ecart-type de l'estimateur de l'aire
    double halfWidth;   // la demi-largeur de l'IC pour l'estimateur de l'aire

    double sum;         // la somme des valeurs
    double sumSquares;  // la somme des carres des valeurs
    uint64_t numGen;      // la taille de l'echantillon (nombre de valeurs generees)

    clock_t start;      // utile pour la mesure du temps requis pour generer un echantillon

public:
    /**
     * Represente le resultat d'un echantillonnage.
     */
    struct Sampling {
        double areaEstimator;                   // aire estimee
        double stdDevEstimator;                 // estimateur de l'ecart-type de l'aire estimee
        ConfidenceInterval confidenceInterval;  // intervalle de confiance a 95%
        uint64_t N;                               // taille de l'echantillon
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

    /**
     * Genere un echantillon d'une taille donnee.
     *
     * @pram N la taille totale de l'echantillon.
     */
    virtual Sampling sampleWithSize(uint64_t N) = 0;

    /**
     * Genere autant de valeurs que necessaire afin d'obtenir une intervalle de confiance a 95% pour l'aire estimee
     * dont la largeur ne depasse pas une certaine valeur donnee.
     *
     * @param maxWidth La taille maximale que doit avoir l'IC.
     * @param step Le nombre de generations qui seront effectuees avant de reverifier la taille de l'IC.
     */
    virtual Sampling sampleWithMaxWidth(double maxWidth, uint64_t step) = 0;

    /**
     * Genere une intervalle de confiance a 95% pour l'aire estimee aussi precise que possible en generant des valeurs
     * durant un laps de temps d'une duree minimum donnee.
     *
     * @param minTime Le temps minimum qui doit etre utilise pour affiner la precision de l'IC.
     * @pram step Le nombre de generations qui seront effectuees avant de reverifier le temps d'execution total.
     */
    virtual Sampling sampleWithMinTime(double minTime, uint64_t step) = 0;

protected:
    /**
     * Initialise les differents champs. Doit etre appelee au debut de chaque etape d'echantillonage.
     */
    void init();
};

#endif // MONTECARLOMETHOD_H
