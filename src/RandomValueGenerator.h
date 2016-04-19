#ifndef RANDOM_VALUE_GENERATOR_H
#define RANDOM_VALUE_GENERATOR_H

#include <random>
#include <vector>
#include "PiecewiseLinearFunction.h"

/**
 * Represente un generateur de realisations de variables aleatoires associees a une fonction affine par morceaux.
 */
class RandomValueGenerator {
protected:
    std::mt19937_64 generator; // generateur mersenne-twister
    std::uniform_real_distribution<double> distribution; // distribution a utiliser pour le mersenne-twister

    PiecewiseLinearFunction func; // la fonction affine par morceaux que l'on utilise
    std::vector<double> F_parts; // parties de la fonction de repartition F

public:
    /**
     * \brief Initialise la fonction par morceaux ainsi que les parties de la fonction de repartition.
     * \param xs Les d'absisses des points constituant la fonction affine par morceaux.
     * \param ys Les ordonnees des points constituant la fonction affine par morceaux.
     */
    RandomValueGenerator(const std::vector<double>& xs, const std::vector<double>& ys);

    /**
     * \brief Initialise la graine du generateur.
     * \param seed La graine a utiliser.
     */
    void setSeed(std::seed_seq seed);

    /**
     * \brief Genere une realisation d'une variable aleatoire associee a la fonction par morceaux.
     */
    virtual double generate() = 0;

protected:
    /**
     * \brief Permet de trouver dans quel intervalle k on tombe en fonction de la probablilité p_k de la tranche liée à
     *        cette intervalle.
     */
    size_t generateK();
};


/**
 * Utilise une approche "bete et mechante" de la methode d'acceptation-rejet afin de generer des realisation de
 * variables aleatoires.
 */
class HitOrMiss : public RandomValueGenerator {
private:
    double a, b; // bornes min et max a prendre en compte pour la generation de l'abcisse du point
    double yMax; // maximum des ordonnees des points constituant "func" (utile pour la génération de l'ordonnee du point)

public:

    /**
     * \brief Initialise les valeurs propres a cet algorithme.
     * \param xs Les d'absisses des points constituant la fonction affine par morceaux.
     * \param ys Les ordonnees des points constituant la fonction affine par morceaux.
     */
    HitOrMiss(const std::vector<double>& xs, const std::vector<double>& ys);

    /**
     * \brief Genere une realisation d'une variable aleatoire associee a la fonction par morceaux.
     */
    double generate();
};


/**
 * Utilise la methode des melanges couplee a une approche geometrique afin de generer des realisation de
 * variables aleatoires.
 */
class Geometric : public RandomValueGenerator {
public:
    /**
     * \brief Initialise les valeurs propres a cet algorithme.
     * \param xs Les d'absisses des points constituant la fonction affine par morceaux.
     * \param ys Les ordonnees des points constituant la fonction affine par morceaux.
     */
    Geometric(const std::vector<double>& xs, const std::vector<double>& ys);

    /**
     * \brief Genere une realisation d'une variable aleatoire associee a la fonction par morceaux.
     */
    double generate();
};


/**
 * Utilise la methode des melanges couplee a la methode des fonctions inverses afin de generer des realisation de
 * variables aleatoires.
 */
class InverseFunctions : public RandomValueGenerator {
public:
    /**
     * \brief Initialise les valeurs propres a cet algorithme.
     * \param xs Les d'abscisses des points constituant la fonction affine par morceaux.
     * \param ys Les ordonnees des points constituant la fonction affine par morceaux.
     */
    InverseFunctions(const std::vector<double>& xs, const std::vector<double>& ys);

    /**
     * \brief Genere une realisation d'une variable aleatoire associee a la fonction par morceaux.
     */
    double generate();
};

#endif // RANDOM_VALUE_GENERATOR_H
