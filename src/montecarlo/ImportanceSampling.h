#ifndef IMPORTANCE_SAMPLING_H
#define IMPORTANCE_SAMPLING_H

#include "MonteCarloMethod.h"
#include "../generators/RandomValueGenerator.h"

/**
 * Represente la methode d'integration par echantillonnage preferentiel.
 */
class ImportanceSampling : public MonteCarloMethod {
private:
    // generateur de variables aleatoires utilisant la methode des melanges couplee a la methode des fonctions inverses
    InverseFunctions generator;

public:
    ImportanceSampling(const Func& g, const std::vector<double>& xs, const std::vector<double>& ys);

    /**
     * @see MontecarloMethod::sampleWithSize.
     */
    Sampling sampleWithSize(uint64_t N);
    /**
     * @see MontecarloMethod::sampleWithMaxWidth.
     */
    Sampling sampleWithMaxWidth(double maxWidth, uint64_t step);
    /**
     * @see MontecarloMethod::sampleWithMinTime.
     */
    Sampling sampleWithMinTime(double maxTime, uint64_t step);

    /**
     * @see MonteCarloMethod::setSeed.
     */
    void setSeed(const std::seed_seq& seed);

private:
    /**
     * Effectue un certain nombre donne de generations afin de mettre a jour les statistiques (somme, somme des
     * carres, moyenne, etc) et de creer un IC pour l'aire estimee.
     *
     * @param step Le nombre de generation qui seront effectuees.
     */
    void sample(uint64_t step);
};

#endif // IMPORTANCE_SAMPLING_H
