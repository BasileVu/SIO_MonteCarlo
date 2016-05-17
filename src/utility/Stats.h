#ifndef STATS_H
#define STATS_H

#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <cstdint>
#include "PiecewiseLinearFunction.h"

/**
 * Represente un intervalle de confiance.
 */
struct ConfidenceInterval {

    double lower; // borne inferieure
    double upper; // borne superieure
    double width; // difference entre les deux bornes
    uint64_t displayPrecision; // precision pour l'affichage

    ConfidenceInterval(double center, double halfWidth, uint64_t precision = 3);

    /**
     * Cree la representation sous forme de chaine de caracteres de l'intervalle de confiance.
     *
     * @return La representation sous forme de chaine de caracteres de l'intervalle de confiance.
     */
    std::string toString() const;
};

/**
 * Permet d'afficher l'intervalle de confiance sur un flux.
 *
 * @param os Le flux de sortie sur lequel ecrire.
 * @pram ci L'intervalle de confiance a utiliser.
 * @return Le flux modifie.
 */
std::ostream& operator<<(std::ostream& os, const ConfidenceInterval& ci);

/**
 * Represente l'ensemble des points d'une fonction affine par morceaux.
 */
struct Points {
    std::vector<double> xs, ys;
};

/**
 * Regroupe differentes fonctions relatives aux statistiques.
 */
class Stats {
public:

    /**
     * Calcule la variance empirique d'un ensemble de donnees.
     *
     * @param values Les donnes a utiliser.
     * @return La variance empirique.
     */
    static double sampleVar(const std::vector<double>& values);

    /**
     * Calcule l'ecart-type empirique d'un ensemble de donnees.
     *
     * @param values Les valeurs a utiliser.
     * @return L'ecart-type empirique.
     */
    static double sampleStdDev(const std::vector<double>& values);

    /**
     * Calcule l'esperance d'une fonction affine par morceaux.
     *
     * @param f La fonction affine par morceaux a utiliser.
     * @return L'esperance.
     */
    static double expectedValue(const PiecewiseLinearFunction& f);

    /**
     * Calcule la moyenne d'un ensemble de donnees.
     *
     * @param values Les valeurs a utiliser.
     * @return La moyenne.
     */
    static double mean(const std::vector<double>& values);

    /**
     * Calcule l'intervalle de confiance d'un ensemble de donnees.
     *
     * @param values Les valeurs a utiliser.
     * @param quantile Le coefficient associe au quantile de l'IC.
     * @return L'intervalle de confiance.
     */
    static ConfidenceInterval confidenceInterval(const std::vector<double>& values, double quantile);

    /**
     * Cree une fonction affine par morceau a partir d'une fonction et d'un nombre de points donnes.
     * Une subdivision reguliere est creee (les largeurs des sous-intervelles sont toutes égales).
     *
     * @param numPoints Le nombre de points qui constitueront la fonction affine par morceaux.
     * @param func La fonction a subdiviser.
     * @param a La borne inferieure de l'intervalle sur lequel on va subdiviser 'func'.
     * @param b La borne superieure de l'intervalle sur lequel on va subdiviser 'func'.
     * @return Les points, un ensemble contenant une liste des abscisses et une des ordonnees.
     */
    static Points createPoints(size_t numPoints, const std::function<double(double)>& func, double a, double b);
};

#endif // STATS_H
