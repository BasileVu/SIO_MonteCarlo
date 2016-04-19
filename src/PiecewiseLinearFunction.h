#ifndef PIECEWISE_LINEAR_FUNCTION_H
#define PIECEWISE_LINEAR_FUNCTION_H

#include <vector>
#include <functional>

/**
 * Regroupe les informations du "morceau" d'une fonction affine par morceaux.
 */
struct Piece {

    // fonction linéaire prenant un double et renvoyant un double
    typedef std::function<double(double)> LinearFunc;

    double x0, x1;    // bornes de l'intervalle definissant le morceau
    double y0, y1;    // ordonnees correspodant aux bornes
    double A_k;       // aire sous la fonction f_k
    LinearFunc f_k;   // fonction associée à ce morceau k.
};

struct PiecewiseLinearFunction {
    std::vector<Piece> pieces;   // "morceaux" de la fonction
    double A = 0;                // aire totale sous la fonction

    PiecewiseLinearFunction(const std::vector<double>& xs, const std::vector<double>& ys);

    /**
     * \brief Recherche dichotomique afin de trouver dans quelle intervalle x se trouve.
     * \param x l'abscisse dont on veut connaître l'intervalle.
     * \return l'indice du morceau dans lequel x se trouve.
     */
    size_t findPart(double x) const;

};

#endif // PIECEWISE_LINEAR_FUNCTION_H
