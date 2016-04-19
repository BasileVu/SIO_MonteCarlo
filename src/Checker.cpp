#include <cstddef>

#include "Checker.h"

bool Checker::check(const std::vector<double>& xs, const std::vector<double>& ys) {

    // verification de la taille des donnees
    if (xs.size() < 2 || ys.size() < 2 || xs.size() != ys.size()) {
        return false;
    }

    return checkXs(xs) && checkYs(ys);
}

bool Checker::checkXs(const std::vector<double>& v) {

    // verification que les abscisses sont strictement croissantes
    for (size_t i = 1; i < v.size(); ++i) {
        if (v[i] <= v[i-1]) {
            return false;
        }
    }

    return true;
}

bool Checker::checkYs(const std::vector<double>& v) {
    bool yNotZero = false;

    // verification que les ordonnees ne sont pas negatives et qu'au moins une ordonnee est plus grande que 0
    for (size_t i = 0; i < v.size(); ++i) {
        if (v[i] < 0) {
            return false;
        }

        if (v[i] > 0) {
            yNotZero = true;
        }
    }

    return yNotZero;
}
