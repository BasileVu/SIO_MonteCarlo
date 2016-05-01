#ifndef CHECKER_H
#define CHECKER_H

#include <vector>

/**
 * Verifie si les donnees sont coherentes afin d'etre utilisees pour generer une fonction affine par morceaux.
 */
class Checker {
public:
    /**
     * \brief Verifie la coherence des donnees.
     * \param xs Les d'abscisses des points constituant la fonction affine par morceaux.
     * \param ys Les ordonnees des points constituant la fonction affine par morceaux.
     * \return Si les donnees sont coherentes.
     */
    static bool check(const std::vector<double>& xs, const std::vector<double>& ys);

private:
    /**
     * \brief Verifie la coherence des donnees pour les abscisses.
     * \param xs Les d'abscisses des points constituant la fonction affine par morceaux.
     * \return Si les donnees sont coherentes.
     */
    static bool checkXs(const std::vector<double>& v);

    /**
     * \brief Verifie la coherence des donnees pour les abscisses.
     * \param ys Les ordonnees des points constituant la fonction affine par morceaux.
     * \return Si les donnees sont coherentes.
     */
    static bool checkYs(const std::vector<double>& v);
};

#endif // CHECKER_H
