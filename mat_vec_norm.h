#ifndef __MAT_VEC_NORM_H_INCLUDED__
#define __MAT_VEC_NORM_H_INCLUDED__

#include "mat_vec.h"

/**
 * @brief Computes the infinite norm of a matrix.
 *        The infinite norm is the maximum absolute row sum.
 *
 * @param M The input matrix.
 * @return The infinite norm of the matrix.
 */
double normInftyMat(const matrix& M);

/**
 * @brief Computes the infinite norm of a banded matrix.
 *        The infinite norm is the maximum absolute row sum.
 *
 * @param M The input banded matrix.
 * @return The infinite norm of the banded matrix.
 */
double normInftyMatBand(const mat_band& M);

/**
 * @brief Computes the Frobenius norm of a matrix.
 *        The Frobenius norm is the square root of the sum of squared matrix elements.
 *
 * @param M The input matrix.
 * @return The Frobenius norm of the matrix.
 */
double normFrobeniusMat(const matrix& M);

/**
 * @brief Computes the Frobenius norm of a banded matrix.
 *        The Frobenius norm is the square root of the sum of squared matrix elements.
 *
 * @param M The input banded matrix.
 * @return The Frobenius norm of the banded matrix.
 */
double normFrobeniusMatBand(const mat_band& M);

/**
 * @brief Computes the infinite norm of a vector.
 *        The infinite norm is the maximum absolute element.
 *
 * @param B The input vector.
 * @return The infinite norm of the vector.
 */
double normInftyVec(const vec& B);

/**
 * @brief Computes the 1-norm of a vector.
 *        The 1-norm is the sum of absolute values of vector elements.
 *
 * @param B The input vector.
 * @return The 1-norm of the vector.
 */
double norm1Vec(const vec& B);

/**
 * @brief Computes the 2-norm of a vector.
 *        The 2-norm is the square root of the sum of squared vector elements.
 *
 * @param B The input vector.
 * @return The 2-norm of the vector.
 */
double norm2Vec(const vec& B);

#endif