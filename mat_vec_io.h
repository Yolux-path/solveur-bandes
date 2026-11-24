// Fourni par l'enseignant.
#ifndef __MAT_VEC_IO_H_INCLUDED__
#define __MAT_VEC_IO_H_INCLUDED__

#include "mat_vec.h"

// Read the coefficients of the matrix from a text file
matrix readMatrix(std::string fileName);
mat_band readMatrixBand(std::string fileName);
mat_band readMatrixBandFromMtx(std::string fileName);

// Write the Band matrix to an ofstream
void writeMatrix(const matrix& M, std::ofstream& output_file);
void writeMatrixBand(const mat_band& M, std::ofstream& output_file);

// Prints on screen the contents of a matrix
void displayMatrix(const matrix& M);
void displayMatrixBand(const mat_band& M);

// Reads the coefficients of a vector from a file
vec readVector(std::string fileName);
vec readVectorFromMtx(std::string fileName);

// Prints on screen the contents of a vector
void displayVector(const vec& B);

// Write the contents of a vector to an ofstream
void writeVector(const vec& B, std::ofstream& output_file);

// Plot Phi = vector of 2D data on a cartesian grid to a txt file for
// gnuplot figures.
void vectorToGnuplot(const vec& Phi, const vec& vx, const vec& vy,
		     std::string fileName, const unsigned int& sample);
#endif
