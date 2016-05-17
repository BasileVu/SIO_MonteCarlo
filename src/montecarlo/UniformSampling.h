#ifndef UNIFORM_SAMPLING_H
#define UNIFORM_SAMPLING_H

#include "MonteCarloMethod.h"

/**
 * Represente la methode d'integration par echantillonnage uniforme.
 */
class UniformSampling : public MonteCarloMethod {
private:
    // generateur mersenne-twister et distribution uniforme
    std::mt19937_64 mtGenerator;
    std::uniform_real_distribution<double> uniformDistr;

    double a, b; // bornes inferieure et superieure de l'intervalle sur lequel on veut evaluer la fonction

public:
    /*
     * Prepare la methode.
     *
     * @param g La fonction dont on veut estimer l'aire.
     * @param a la borne inferieure de l'intervalle sur lequel on veut evaluer g.
     * @param b la borne superieure de l'intervalle sur lequel on veut evaluer g.
     */
    UniformSampling(const Func& g, double a, double b);


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
     * carres, moyenne, etc) afin de creer un IC pour l'aire estimee.
     *
     * @param step Le nombre de generation qui seront effectuees.
     */
    void sample(uint64_t step);

    /**
     * Empaquete toutes les valeurs associees a l'echantillon (aire estimee, IC, etc: voir MonteCarloMethod::Sampling).
     *
     * @return L'echantillon cree.
     */
    Sampling createSampling(double timeElapsed) const;
};

#endif // UNIFORMS_AMPLING_H
