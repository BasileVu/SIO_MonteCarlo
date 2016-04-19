#include <algorithm>

#include "Checker.h"
#include "RandomValueGenerator.h"

RandomValueGenerator::RandomValueGenerator(const std::vector<double>& xs, const std::vector<double>& ys)
            : func(xs, ys), distribution(std::uniform_real_distribution<double>(0, 1)) {

    // verification de la coherence des donnees
    if (!Checker::check(xs, ys)) {
        throw std::invalid_argument("Erreur: Les donnees ne sont pas coherentes.");
    }

    // creation des pk
    std::vector<double> pks;
    pks.reserve(xs.size()-1);
    for (size_t i = 0; i < xs.size()-1; ++i) {
        pks.push_back(func.pieces[i].A_k/func.A);
    }

    // préparation des parties de F
    F_parts.resize(xs.size());
    F_parts[0] = 0; // F_0 : premiere partie de la fonction de repartition -> 0 avant xs[0]

    for (size_t i = 1; i < F_parts.size(); ++i) {
        F_parts[i] = F_parts[i-1] + pks[i-1];
    }
}

void RandomValueGenerator::setSeed(std::seed_seq seed) {
    generator.seed(seed);
}

size_t RandomValueGenerator::generateK() {
    size_t j = 1;
    double U = distribution(generator);

    // on cherche l'indice de l'intervalle dans lequel on est tombe
    while (true) {
        if (U <= F_parts[j]) {
            return j-1;
        }
        ++j;
    }
}

HitOrMiss::HitOrMiss(const std::vector<double>& xs, const std::vector<double>& ys)
        : RandomValueGenerator(xs, ys) {
    a = xs.front(), b = xs.back();
    yMax = *std::max_element(ys.begin(), ys.end());
}

Geometric::Geometric(const std::vector<double>& xs, const std::vector<double>& ys)
        : RandomValueGenerator(xs, ys) {}

InverseFunctions::InverseFunctions(const std::vector<double>& xs, const std::vector<double>& ys)
        : RandomValueGenerator(xs, ys) {}


double HitOrMiss::generate() {

    double X, Y; // coordonnees du point (X,Y) qui sera genere
    size_t sliceIndex;  // indice de la "tranche" dans laquelle X se trouvera

    do {
        // generation du point (X,Y)
        X = distribution(generator) * (b - a) + a;
        Y = distribution(generator) * yMax;

        // on cherche le morceau lie à l'intervalle dans laquelle X se trouve
        sliceIndex = func.findPart(X);

        // rejet si Y est > que f(X), avec f_k la fonction affine associee a l'intervalle k
    } while (Y > func.pieces[sliceIndex].f_k(X));

    return X;
}


double Geometric::generate() {

    // On commence par selectionner une intervalle en fonction des p_k des "tranches" de la fonction.
    // K represente l'indice de l'intervalle selectionne.
    size_t K = generateK();


    // Ensuite, on genere une realisation d'une variable de densite f_K en acceptant a tous les coups X.

    // morceau associe a l'indice K
    const Piece& piece = func.pieces[K];

    // valeurs relatives au morceau de fonction K
    double x0 = piece.x0, x1 = piece.x1;
    double yMax = std::max(piece.y0, piece.y1);

    // generation du point (X,Y)
    double X = distribution(generator) * (x1 - x0) + x0;
    double Y = distribution(generator) * yMax;

    // Si Y est sous f_K, ok, on retourne X. Sinon, on applique une symetrie à X et on le retourne.
    if (Y <= piece.f_k(X)) {
        return X;
    } else {
        return x0 + (x1 - X);
    }
}

double InverseFunctions::generate() {

    // On commence par selectionner une intervalle en fonction des p_k des "tranches" de la fonction.
    // K represente l'indice de l'intervalle selectionne.
    size_t K = generateK();

    // Ensuite, on applique la methode des fonctions inverses.

    // morceau associe a l'indice K
    const Piece& piece = func.pieces[K];

    // valeurs relatives au morceau de fonction K
    double x0 = piece.x0, x1 = piece.x1;
    double y0 = piece.y0, y1 = piece.y1;

    double U = distribution(generator);

    // Si y0 = y1, alors on est dans le cas d'une uniforme. Sinon, on inverse la fonction de repartition associee
    // a la fonction f_K.
    if (y0 == y1) {
        return x0 + U*(x1 - x0);
    } else {
        double m = (y1 - y0)/(x1 - x0);
        return x0 + (sqrt( (y1*y1 - y0*y0) * U + y0*y0 ) - y0) / m;
    }
}
