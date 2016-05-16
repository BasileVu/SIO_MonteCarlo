#ifndef CONTROL_VARIABLE_H
#define CONTROL_VARIABLE_H

#include "MonteCarloMethod.h"

/**
 * Represente la methode d'integration par echantillonage uniforme avec variable de controle.
 */
class ControlVariable : public MonteCarloMethod {
private:
    std::mt19937_64 mtGenerator;
    std::uniform_real_distribution<double> uniformDistr;

    PiecewiseLinearFunction h; // variable de controle (fonction affine par morceaux)

    double a, b; // bornes inferieure et superieure
    double mu;   // esperance de la fonction par morceaux

    double c; // coefficient c tq V = Y + c(Z - mu), avec Y = g(X) et Z(X) la variable de controle

public:
    /**
     * Prepare la methode.
     *
     * @param g La fonction dont on veut connaitre l'aire.
     * @param a la borne inferieure de l'intervalle sur laquelle on veut evaluer g.
     * @param b la borne superieure de l'intervalle sur laquelle on veut evaluer g.
     * @param xs Les abscisses des points de la variable de controle (sous forme de fonciton par morceaux).
     * @param ys Les ordonnees des points de la variable de controle (sous forme de fonciton par morceaux).
     */
    ControlVariable(const Func& g, double a, double b, const std::vector<double>& xs, const std::vector<double>& ys);

    /**
     * Genere un echantillon d'une taille donnee.
     *
     * @param M La taille de l'echantillon "de base" (pour trouver le 'c').
     * @pram N la taille totale de l'echantillon.
     */
    Sampling sampleWithSize(size_t M, size_t N);

    /**
     * Genere autant de valeurs que necessaire afin d'obtenir une intevralle de confiance a 95% pour l'aire estimee
     * dont la largeur ne depasse pas une certaine valeur donnee.
     *
     * @param M La taille de l'echantillon "de base" (pour trouver le 'c').
     * @param
     */
    Sampling sampleWithMaxWidth(size_t M, double maxWidth, size_t step);

    Sampling sampleWithMinTime(size_t M, double maxTime, size_t step);

    void setSeed(const std::seed_seq& seed);

private:
    void computeConstant(size_t M);

    void sample(size_t step);

    Sampling createSampling(double timeElapsed) const;
};

#endif // CONTROL_VARIABLE_H
