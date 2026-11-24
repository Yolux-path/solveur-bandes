#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "linalg_band.h"

using uint = unsigned int;

namespace {

// Bornes [debut, fin) des lignes non nulles de la colonne j
uint rowBegin(const mat_band& A, uint j) { return j > A.ku ? j - A.ku : 0; }
uint rowEnd(const mat_band& A, uint j) { return std::min(A.n, j + A.kl + 1); }

double& at(mat_band& A, uint i, uint j) { return A.col[j][A.ku + i - j]; }
double at(const mat_band& A, uint i, uint j) { return A.col[j][A.ku + i - j]; }

}  // namespace

double traceBand(const mat_band& M) {
    double result = 0;
    for (uint j = 0; j < M.n; j++) result += at(M, j, j);
    return result;
}

vec matBandVec(const mat_band& M, const vec& X) {
    vec Y(M.n, 0.);
    for (uint j = 0; j < M.n; j++)
        for (uint i = rowBegin(M, j); i < rowEnd(M, j); i++) Y[i] += at(M, i, j) * X[j];
    return Y;
}

vec residualBand(const vec& B, const mat_band& M, const vec& X) {
    vec R = matBandVec(M, X);
    for (uint i = 0; i < M.n; i++) R[i] = B[i] - R[i];
    return R;
}

bool isSymmetric(const mat_band& M) {
    for (uint j = 0; j < M.n; j++)
        for (uint i = rowBegin(M, j); i < rowEnd(M, j); i++) {
            // A(j,i) est nul s'il est hors de la bande
            double t = (j + M.ku >= i && i + M.kl >= j) ? at(M, j, i) : 0.;
            double a = at(M, i, j);
            if (std::abs(a - t) > 1e-12 * std::max(std::abs(a), std::abs(t))) return false;
        }
    return true;
}

int CholeskiBand(mat_band& A) {
    if (!isSymmetric(A)) return -1;

    for (uint j = 0; j < A.n; j++) {
        for (uint k = (j > A.kl ? j - A.kl : 0); k < j; k++) {
            const double ljk = at(A, j, k);
            for (uint i = j; i < rowEnd(A, k); i++) at(A, i, j) -= at(A, i, k) * ljk;
        }
        if (!(at(A, j, j) > 0.)) return -1;
        at(A, j, j) = std::sqrt(at(A, j, j));
        for (uint i = j + 1; i < rowEnd(A, j); i++) at(A, i, j) /= at(A, j, j);
        for (uint i = rowBegin(A, j); i < j; i++) at(A, i, j) = 0.;
    }
    return 0;
}

vec CholeskiBand_solve(const mat_band& L, const vec& B) {
    vec X = B;
    for (uint j = 0; j < L.n; j++) {
        X[j] /= at(L, j, j);
        for (uint i = j + 1; i < rowEnd(L, j); i++) X[i] -= at(L, i, j) * X[j];
    }
    for (uint j = L.n; j-- > 0;) {
        for (uint i = j + 1; i < rowEnd(L, j); i++) X[j] -= at(L, i, j) * X[i];
        X[j] /= at(L, j, j);
    }
    return X;
}

double CholeskiBand_det(const mat_band& L) {
    double det = 1;
    for (uint i = 0; i < L.n; i++) det *= at(L, i, i) * at(L, i, i);
    return det;
}

int LUBand(mat_band& A, std::vector<uint>& piv) {
    // Place pour le remplissage : kl diagonales supérieures de plus
    for (auto& c : A.col) c.insert(c.begin(), A.kl, 0.);
    A.ku += A.kl;
    piv.assign(A.n, 0);

    for (uint k = 0; k < A.n; k++) {
        uint p = k;
        for (uint i = k + 1; i < rowEnd(A, k); i++)
            if (std::abs(at(A, i, k)) > std::abs(at(A, p, k))) p = i;
        piv[k] = p;
        if (at(A, p, k) == 0.) return -1;

        const uint jEnd = std::min(A.n, k + A.ku + 1);
        if (p != k)
            for (uint j = k; j < jEnd; j++) std::swap(at(A, k, j), at(A, p, j));

        for (uint i = k + 1; i < rowEnd(A, k); i++) at(A, i, k) /= at(A, k, k);
        for (uint j = k + 1; j < jEnd; j++) {
            const double akj = at(A, k, j);
            if (akj == 0.) continue;
            for (uint i = k + 1; i < rowEnd(A, k); i++) at(A, i, j) -= at(A, i, k) * akj;
        }
    }
    return 0;
}

vec LUBand_solve(const mat_band& A, const std::vector<uint>& piv, const vec& B) {
    vec X = B;
    for (uint k = 0; k < A.n; k++) {
        std::swap(X[k], X[piv[k]]);
        for (uint i = k + 1; i < rowEnd(A, k); i++) X[i] -= at(A, i, k) * X[k];
    }
    for (uint j = A.n; j-- > 0;) {
        X[j] /= at(A, j, j);
        for (uint i = rowBegin(A, j); i < j; i++) X[i] -= at(A, i, j) * X[j];
    }
    return X;
}

double LUBand_det(const mat_band& A, const std::vector<uint>& piv) {
    double det = 1;
    for (uint i = 0; i < A.n; i++) det *= (piv[i] != i ? -1 : 1) * at(A, i, i);
    return det;
}

double logAbsDiag(const mat_band& A) {
    double s = 0;
    for (uint i = 0; i < A.n; i++) s += std::log(std::abs(at(A, i, i)));
    return s;
}

int SOR(const mat_band& A, const vec& b, vec& x, double omega, int max_iter, double tolerance) {
    for (uint i = 0; i < A.n; i++)
        if (at(A, i, i) == 0.) {
            std::cerr << "SOR: zero diagonal entry in row " << i << std::endl;
            return max_iter;
        }

    for (int iter = 1; iter <= max_iter; iter++) {
        double change = 0.;
        for (uint i = 0; i < A.n; i++) {
            double sum = b[i];
            for (uint j = (i > A.kl ? i - A.kl : 0); j < std::min(A.n, i + A.ku + 1); j++)
                if (j != i) sum -= at(A, i, j) * x[j];
            const double xi = (1.0 - omega) * x[i] + omega * sum / at(A, i, i);
            change += std::abs(xi - x[i]);
            x[i] = xi;
        }
        if (!std::isfinite(change)) break;
        if (change < tolerance) return iter;
    }
    std::cerr << "SOR did not converge." << std::endl;
    return max_iter;
}
