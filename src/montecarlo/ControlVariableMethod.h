#ifndef CONTROL_VARIABLE_H
#define CONTROL_VARIABLE_H

#include "MonteCarloMethod.h"

/**
 * Represente la methode d'integration par echantillonage uniforme avec variable de controle.
 */
class ControlVariable : public MonteCarloMethod {
private:
    // generateur mersenne-twister et distribution uniforme
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
     * Genere autant de valeurs que necessaire afin d'obtenir une intervalle de confiance a 95% pour l'aire estimee
     * dont la largeur ne depasse pas une certaine valeur donnee.
     *
     * @param M La taille de l'echantillon "de base" (pour trouver le 'c').
     * @param maxWidth La taille maximale que doit avoir l'IC.
     * @pram step Le nombre de generations qui seront effectuees avant de reverifier la taille de l'IC.
     */
    Sampling sampleWithMaxWidth(size_t M, double maxWidth, size_t step);

    /**
     * Genere une intervalle de confiance a 95% pour l'aire estimee aussi precise que possible en generant des valeurs
     * durant un laps de temps d'une duree minimum donnee.
     *
     * @param M La taille de l'echantillon "de base" (pour trouver le 'c').
     * @param minTime Le temps minimum qui doit etre utilise pour affiner la precision de l'IC.
     * @pram step Le nombre de generations qui seront effectuees avant de reverifier le temps d'execution total.
     */
    Sampling sampleWithMinTime(size_t M, double minTime, size_t step);

    /**
     * Initialise la graine du generateur.
     *
     * @param seed La graine a utiliser.
     */
    void setSeed(const std::seed_seq& seed);

private:
    /**
     * Calcule la constante 'c'.
     *
     * @param M La taille de l'echantillon a generer afin de toruver la constante.
     */
    void computeConstant(size_t M);

    /**
     * Utilise la constante generee afin de generer l'IC. Effectue un certain nombre donne de generations afin de mettre
     * a jour les statistiques (somme et somme des carres) et de creer un IC.
     *
     * @param step Le nomre de generation qui seront effectuees.
     */
    void sample(size_t step);

    /**
     * Cree l'echantillon.
     *
     * @return L'echantillon cree.
     */
    Sampling createSampling(double timeElapsed) const;
};

#endif // CONTROL_VARIABLE_H
