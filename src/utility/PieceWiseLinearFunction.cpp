#include "PiecewiseLinearFunction.h"

PiecewiseLinearFunction::PiecewiseLinearFunction(const std::vector<double>& xs, const std::vector<double>& ys) {

    pieces.reserve(xs.size());
    for (size_t i = 0; i < xs.size(); ++i) {

        // creation d'un morceau de la fonction
        if (i < xs.size() - 1) {

            // ajout preliminaire des abscisses et ordonnees
            Piece s {xs[i], xs[i+1], ys[i], ys[i+1]};

            // aire sous le morceau de fonction
            s.A_k = (ys[i+1] + ys[i]) * (xs[i+1] - xs[i]) / 2;

            // morceau de fonction
            s.f_k = [&xs, &ys, i](double x) {
                double m = (ys[i + 1] - ys[i]) / (xs[i+1] - xs[i]);
                return m * (x - xs[i]) + ys[i];
            };

            // ajout du morceau
            pieces.push_back(s);
            A += s.A_k;
        }
    }
}

size_t PiecewiseLinearFunction::findPart(double x) const {

    size_t first = 0, last = pieces.size() - 1; // indices des tranches à prendre en compte

    while (true) {

        // il ne reste qu'une tranche -> x est dedans
        if (first == last) {
            return first;
        }

        // on regarde si x est dans la première ou deuxième moitié de l'intervalle de recherche
        size_t mid = (last-first)/2 + first;
        if (x < pieces[mid].x1) {   // première moitié
            last = mid;
        } else {                    // deuxième moitié
            first = mid+1;
        }
    }
}