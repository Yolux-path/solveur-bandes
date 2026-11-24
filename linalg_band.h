#ifndef __LINALG_BAND_H_INCLUDED__
#define __LINALG_BAND_H_INCLUDED__

#include "mat_vec.h"

// Stockage bande : A(i,j) = col[j][ku + i - j] pour -ku <= i - j <= kl.

/** @brief Trace de la matrice. */
double traceBand(const mat_band& M);

/** @brief Produit Y = M * X. */
vec matBandVec(const mat_band& M, const vec& X);

/** @brief Résidu B - M * X. */
vec residualBand(const vec& B, const mat_band& M, const vec& X);

/** @brief Vrai si M est symétrique (à 1e-12 près en relatif). */
bool isSymmetric(const mat_band& M);

/**
 * @brief Factorisation de Cholesky A = L * L^T, en place (L dans la partie inférieure).
 * @return 0, ou -1 si A n'est pas symétrique définie positive.
 */
int CholeskiBand(mat_band& A);

/** @brief Résout L * L^T * X = B. */
vec CholeskiBand_solve(const mat_band& L, const vec& B);

/** @brief Déterminant à partir de la factorisation de Cholesky. */
double CholeskiBand_det(const mat_band& L);

/**
 * @brief Factorisation LU en place avec pivotage partiel.
 * La bande supérieure est élargie à ku + kl pour le remplissage dû aux échanges.
 * @param piv piv[k] = ligne échangée avec la ligne k.
 * @return 0, ou -1 si la matrice est singulière.
 */
int LUBand(mat_band& A, std::vector<unsigned int>& piv);

/** @brief Résout A * X = B à partir de la factorisation LU et des pivots. */
vec LUBand_solve(const mat_band& A, const std::vector<unsigned int>& piv, const vec& B);

/** @brief Déterminant à partir de la factorisation LU et des pivots. */
double LUBand_det(const mat_band& A, const std::vector<unsigned int>& piv);

/** @brief log|det| à partir d'une factorisation (somme des logs, sans débordement). */
double logAbsDiag(const mat_band& A);

/**
 * @brief Méthode SOR ; x contient l'itéré initial puis la solution.
 * @return Nombre d'itérations, ou max_iter si pas de convergence.
 */
int SOR(const mat_band& A, const vec& b, vec& x, double omega, int max_iter, double tolerance);

#endif
