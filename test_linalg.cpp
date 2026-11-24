// Compare les routines bande à un calcul dense, sur des matrices aléatoires.
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>

#include "linalg_band.h"
#include "mat_vec_norm.h"

#define CHECK(c) do { if (!(c)) { std::cerr << "ECHEC ligne " << __LINE__ << " : " #c << std::endl; std::exit(1); } } while (0)

using Dense = std::vector<vec>;

static std::mt19937 gen(42);
static double rnd() { return std::uniform_real_distribution<double>(-1, 1)(gen); }

static mat_band make(unsigned n, unsigned kl, unsigned ku) {
    mat_band A{n, kl, ku, std::vector<vec>(n, vec(kl + ku + 1, 0.))};
    return A;
}
static bool inBand(const mat_band& A, int i, int j) { return i - j <= int(A.kl) && j - i <= int(A.ku); }
static double& el(mat_band& A, int i, int j) { return A.col[j][A.ku + i - j]; }
static Dense dense(mat_band A) {
    Dense D(A.n, vec(A.n, 0.));
    for (int i = 0; i < int(A.n); i++)
        for (int j = 0; j < int(A.n); j++)
            if (inBand(A, i, j)) D[i][j] = el(A, i, j);
    return D;
}
static double dist(const vec& a, const vec& b) {
    double m = 0;
    for (size_t i = 0; i < a.size(); i++) m = std::max(m, std::abs(a[i] - b[i]));
    return m;
}
static double detDense(Dense D) {
    double det = 1;
    int n = D.size();
    for (int k = 0; k < n; k++) {
        int p = k;
        for (int i = k; i < n; i++) if (std::abs(D[i][k]) > std::abs(D[p][k])) p = i;
        if (p != k) { std::swap(D[p], D[k]); det = -det; }
        det *= D[k][k];
        for (int i = k + 1; i < n; i++)
            for (int j = n - 1; j >= k; j--) D[i][j] -= D[i][k] / D[k][k] * D[k][j];
    }
    return det;
}

int main() {
    for (unsigned n : {1u, 2u, 5u, 20u, 60u})
        for (unsigned kl : {0u, 1u, 3u})
            for (unsigned ku : {0u, 2u, 4u}) {
                if (kl >= n || ku >= n) continue;
                mat_band A = make(n, kl, ku);
                vec x(n);
                for (int i = 0; i < int(n); i++) {
                    x[i] = rnd();
                    for (int j = 0; j < int(n); j++)
                        // diagonale nulle une fois sur trois : impose le pivotage
                        if (inBand(A, i, j)) el(A, i, j) = (i == j && kl > 0 && ku > 0 && n > 2 && i % 3 == 0) ? 0. : rnd();
                }
                Dense D = dense(A);
                vec b(n, 0.);
                double tr = 0;
                for (unsigned i = 0; i < n; i++) {
                    tr += D[i][i];
                    for (unsigned j = 0; j < n; j++) b[i] += D[i][j] * x[j];
                }
                CHECK(dist(matBandVec(A, x), b) < 1e-13);
                CHECK(std::abs(traceBand(A) - tr) < 1e-13);

                mat_band F = A;
                std::vector<unsigned> piv;
                CHECK(LUBand(F, piv) == 0);
                vec y = LUBand_solve(F, piv, b);
                CHECK(normInftyVec(residualBand(b, A, y)) < 1e-11 * (1 + normInftyVec(y)));
                double det = LUBand_det(F, piv), ref = detDense(D);
                CHECK(std::abs(det - ref) <= 1e-9 * std::abs(ref));
                CHECK(std::abs(logAbsDiag(F) - std::log(std::abs(ref))) < 1e-8 * (1 + std::abs(std::log(std::abs(ref)))));

                // Matrice symétrique à diagonale strictement dominante : SDP
                mat_band S = make(n, kl, kl);
                for (int i = 0; i < int(n); i++) {
                    el(S, i, i) = 2. * kl + 1.;
                    for (int j = std::max(0, i - int(kl)); j < i; j++) el(S, i, j) = el(S, j, i) = rnd();
                }
                CHECK(isSymmetric(S));
                vec bs = matBandVec(S, x);
                mat_band L = S;
                CHECK(CholeskiBand(L) == 0);
                CHECK(dist(CholeskiBand_solve(L, bs), x) < 1e-12);
                Dense Ld = dense(L), Sd = dense(S);
                for (unsigned i = 0; i < n; i++)
                    for (unsigned j = 0; j < n; j++) {
                        double s = 0;
                        for (unsigned k = 0; k < n; k++) s += Ld[i][k] * Ld[j][k];
                        CHECK(std::abs(s - Sd[i][j]) < 1e-12 * (2. * kl + 1.));
                    }
                CHECK(std::abs(CholeskiBand_det(L) / detDense(Sd) - 1) < 1e-10);
                for (double w : {0.8, 1.0, 1.3}) {
                    vec z(n, 0.);
                    CHECK(SOR(S, bs, z, w, 10000, 1e-13) < 10000);
                    CHECK(dist(z, x) < 1e-11);
                }

                // Non symétrique à diagonale dominante : Gauss-Seidel converge
                mat_band N = A;
                for (int i = 0; i < int(n); i++) el(N, i, i) = kl + ku + 1.;
                vec bn = matBandVec(N, x), z(n, 0.);
                CHECK(SOR(N, bn, z, 1.0, 10000, 1e-13) < 10000);
                CHECK(dist(z, x) < 1e-11);
            }

    // Cas limites
    std::vector<unsigned> piv;
    mat_band P = make(2, 1, 1);
    el(P, 0, 1) = el(P, 1, 0) = 1.;
    CHECK(isSymmetric(P));
    mat_band C = P;
    CHECK(CholeskiBand(C) == -1);  // symétrique mais pas définie positive
    CHECK(LUBand(P, piv) == 0);
    CHECK(dist(LUBand_solve(P, piv, {3., 7.}), {7., 3.}) < 1e-15);
    CHECK(LUBand_det(P, piv) == -1.);

    mat_band Z = make(3, 1, 1);
    CHECK(LUBand(Z, piv) == -1);  // singulière
    vec z(3, 0.);
    CHECK(SOR(Z, {1., 1., 1.}, z, 1., 10, 1e-6) == 10);  // diagonale nulle

    mat_band U = make(3, 0, 1);
    el(U, 0, 1) = 1.;
    CHECK(!isSymmetric(U));
    CHECK(std::isnan(normInftyVec({1., std::nan(""), 2.})));
    CHECK(norm1Vec({-1., 2., -3.}) == 6.);

    std::cout << "Tous les tests passent." << std::endl;
}
