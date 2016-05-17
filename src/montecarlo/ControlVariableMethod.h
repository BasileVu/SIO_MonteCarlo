#ifndef CONTROL_VARIABLE_H
#define CONTROL_VARIABLE_H

#include "MonteCarloMethod.h"

/**
 * Represente la methode d'integration par echantillonnage uniforme avec variable de controle.
 */
class ControlVariable : public MonteCarloMethod {
private:
    // generateur mersenne-twister et distribution uniforme
    std::mt19937_64 mtGenerator;
    std::uniform_real_distribution<double> uniformDistr;

    PiecewiseLinearFunction h; // variable de controle (fonction affine par morceaux)

    double a, b; // bornes inferieure et superieure de l'intervalle sur laquelle on veut evaluer la fonction
    double mu;   // esperance de la fonction par morceaux

    double c;     // coefficient c tq V = Y + c(Z - mu), avec Y = g(X) et Z(X) la variable de controle
    uint64_t M = 0; // taille de l'echantillon pour determiner 'c'

public:
    /**
     * Prepare la methode.
     *
     * @param g La fonction dont on veut estimer l'aire.
     * @param a la borne inferieure de l'intervalle sur laquelle on veut evaluer g.
     * @param b la borne superieure de l'intervalle sur laquelle on veut evaluer g.
     * @param xs Les abscisses des points de la variable de controle (sous forme de fonciton par morceaux).
     * @param ys Les ordonnees des points de la variable de controle (sous forme de fonciton par morceaux).
     */
    ControlVariable(const Func& g, double a, double b, const std::vector<double>& xs, const std::vector<double>& ys);

    /**
     * Fixe la valeur de
     */
    void setSamplingSize(uint64_t M);


    /**
     * @see MontecarloMethod::sampleWithSize.
     *
     * Utilisable uniquement apres que 'setSamplingSize' ait ete appelee au moins une fois apres la creation de l'objet.
     */
    Sampling sampleWithSize(uint64_t N);
    /**
     * @see MontecarloMethod::sampleWithMaxWidth.
     *
     * Utilisable uniquement apres que 'setSamplingSize' ait ete appelee au moins une fois apres la creation de l'objet.
     */
    Sampling sampleWithMaxWidth(double maxWidth, uint64_t step);
    /**
     * @see MontecarloMethod::sampleWithMinTime.
     *
     * Utilisable uniquement apres que 'setSamplingSize' ait ete appelee au moins une fois apres la creation de l'objet.
     */
    Sampling sampleWithMinTime(double maxTime, uint64_t step);

    /**
     * @see MonteCarloMethod::setSeed.
     */
    void setSeed(const std::seed_seq& seed);

private:
    /**
     * Calcule la constante 'c'.
     */
    void computeConstant();

    /*
     * Effectue un certain nombre donne de generations afin de mettre a jour les statistiques (somme, somme des
     * carres, moyenne, etc) afin de creer un IC pour l'aire estimee.
     *
     * Utilisable uniquement apres un appel a 'computeConstant'.
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

#endif // CONTROL_VARIABLE_H
