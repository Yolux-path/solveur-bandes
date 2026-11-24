#include <string>
#include <iostream> // The std library for i/o
#include <fstream>  // For file streams
#include <sstream>  // For string streams
#include <iomanip>  // Manipulate i/o, e.g. formatting
#include <cmath>
#include <ctime>
#include <limits>
#include <vector>

#include "mat_vec.h" // My struct for matrices
#include "linalg_band.h"
#include "mat_vec_norm.h" //My functions who computes normes

constexpr double EPSILON = 1e-20;


//Write the line (a): n kl ku nnz(j<i) nnz(j=i) nnz(j>i) nnz(total) symetrique(O/N)
std::vector<uint> ComputeAInfos(const mat_band& M) {
    uint countnnz = 0;
    uint countnnztotal = 0;
    std::vector<uint> AInfos(8);

    AInfos[0] = M.n;
    AInfos[1] = M.kl;
    AInfos[2] = M.ku;

    for (uint j = 0; j < M.n; j++) {
        for (uint k = M.ku + 1; k < M.ku + M.kl + 1; k++) {
            if (j + k >= M.ku && j + k < M.n + M.ku) {
                if (std::abs(M.col[j][k]) > EPSILON) {
                    countnnz++;
                }
            }
        }
    }

    AInfos[3] = countnnz;
    countnnztotal += countnnz;

    countnnz = 0;
    for (uint j = 0; j < M.n; j++) {
        if (std::abs(M.col[j][M.ku]) > EPSILON) {
            countnnz++;
        }
    }
    AInfos[4] = countnnz;
    countnnztotal += countnnz;

    countnnz = 0;
    for (uint j = 0; j < M.n; j++) {
        for (uint k = 0; k < M.ku; k++) {
            if (j + k >= M.ku && j + k < M.n + M.ku) {
                if (std::abs(M.col[j][k]) > EPSILON) {
                    countnnz++;
                }
            }
        }
    }
    AInfos[5] = countnnz;
    countnnztotal += countnnz;

    AInfos[6] = countnnztotal;

    bool symetrique = isSymmetric(M);

    if (symetrique) {
        AInfos[7] = 1;
    } else {
        AInfos[7] = 0;
    }

    return AInfos;
}

// Function to write A table
void writeMatrixBandAtable(std::ofstream& output_file, const std::vector<uint>& Infos) {
    output_file << "# a" << std::endl;
    output_file << std::setw(15) << "n : " << Infos[0] << std::endl;
    output_file << std::setw(15) << "kl : " << Infos[1] << std::endl;
    output_file << std::setw(15) << "ku : " << Infos[2] << std::endl;
    output_file << std::setw(15) << "nnz (j < i ) : " << Infos[3] << std::endl;
    output_file << std::setw(15) << "nnz (j = i ) : " << Infos[4] << std::endl;
    output_file << std::setw(15) << "nnz (j > i ) : " << Infos[5] << std::endl;
    output_file << std::setw(15) << "nnz (total ) : " << Infos[6] << std::endl;
    output_file << std::setw(15) << "symetrique (O/N) : " << (Infos[7] == 1 ? "Symetrique" : "Non Symetrique") << std::endl << std::endl;
}

//Write the line (b): <  =  > Dominante(O/N) DiagonaleDominante(O/N) DefiniePositive(O/N)
std::vector<uint> ComputeBInfos(const mat_band& M) {
    std::vector<uint> BInfos(4, 0);
    std::vector<double> summ(M.n, 0.0);

    for (uint j = 0; j < M.n; j++) {
        for (uint k = 0; k < M.ku + M.kl + 1; k++) {
            if (j + k >= M.ku && j + k < M.n + M.ku && k != M.ku) {
                summ[k + j - M.ku] += std::abs(M.col[j][k]);
            }
        }
    }

    for (int j = 0; j < (int)M.n; j++) {
        const double diag = std::abs(M.col[j][M.ku]);
        if (std::abs(summ[j] - diag) <= 1e-12 * diag) {
            BInfos[1]++;
        } else if (summ[j] > diag) {
            BInfos[0]++;
        } else {
            BInfos[2]++;
        }
    }

    // 0 non; 1 oui; 2 strictement; 3 fortement (au moins une ligne stricte)
    if (BInfos[0] > 0) {
        BInfos[3] = 0;
    } else if (BInfos[2] == M.n) {
        BInfos[3] = 2;
    } else if (BInfos[2] > 0) {
        BInfos[3] = 3;
    } else {
        BInfos[3] = 1;
    }

    return BInfos;
}

// Function to write B table
void writeMatrixBandBtable(std::ofstream& output_file, const std::vector<uint>& Infos) {
    output_file << "# b" << std::endl;
    output_file << "Card(i : abs(Aij) < somme(abs(Aij)) j≠i) : " << Infos[0] << std::endl;
    output_file << "Card(i : abs(Aij) = somme(abs(Aij)) j≠i) : " << Infos[1] << std::endl;
    output_file << "Card(i : abs(Aij) > somme(abs(Aij)) j≠i) : " << Infos[2] << std::endl;
    output_file << "Dominante (non; oui; fortement; strictement) : " << std::endl;
    switch (Infos[3]) {
        case 0: output_file << "Non Dominante"; break;
        case 1: output_file << "Dominante"; break;
        case 2: output_file << "Strictement Dominante"; break;
        case 3: output_file << "Fortement Dominante"; break;
    }
    output_file << std::endl << std::endl;
}


//Write the line (c): min(j<i) max(j<i) min(j=i) max(j=i) min(j>i) max(j>i)
std::vector<double> ComputeCInfos(const mat_band& M) {
    std::vector<double> CInfo(6);

    // Find the smallest value greater than EPSILON (j < i)
    CInfo[0] = std::numeric_limits<double>::max();
    for (uint j = 0; j < M.n; j++) {
        for (uint k = M.ku + 1; k < M.ku + M.kl + 1; k++) {
            if (j + k >= M.ku && j + k < M.n + M.ku) {
                if (std::abs(M.col[j][k]) >= EPSILON && std::abs(M.col[j][k]) < CInfo[0]) {
                    CInfo[0] = std::abs(M.col[j][k]);
                }
            }
        }
    }

    // Find the largest value (j < i)
    CInfo[1] = 0;
    for (uint j = 0; j < M.n; j++) {
        for (uint k = M.ku + 1; k < M.ku + M.kl + 1; k++) {
            if (j + k >= M.ku && j + k < M.n + M.ku) {
                if (std::abs(M.col[j][k]) > CInfo[1]) {
                    CInfo[1] = std::abs(M.col[j][k]);
                }
            }
        }
    }

    // Find the smallest value greater than EPSILON (j = i)
    CInfo[2] = std::numeric_limits<double>::max();
    for (uint j = 0; j < M.n; j++) {
        if (std::abs(M.col[j][M.ku]) >= EPSILON && std::abs(M.col[j][M.ku]) < CInfo[2]) {
            CInfo[2] = std::abs(M.col[j][M.ku]);
        }
    }

    // Find the largest value (j = i)
    CInfo[3] = 0;
    for (uint j = 0; j < M.n; j++) {
        if (std::abs(M.col[j][M.ku]) > CInfo[3]) {
            CInfo[3] = std::abs(M.col[j][M.ku]);
        }
    }

    // Find the smallest value greater than EPSILON (j > i)
    CInfo[4] = std::numeric_limits<double>::max();
    for (uint j = 0; j < M.n; j++) {
        for (uint k = 0; k < M.ku; k++) {
            if (j + k >= M.ku && j + k < M.n + M.ku) {
                if (std::abs(M.col[j][k]) >= EPSILON && std::abs(M.col[j][k]) < CInfo[4]) {
                    CInfo[4] = std::abs(M.col[j][k]);
                }
            }
        }
    }

    // Find the largest value (j > i)
    CInfo[5] = 0;
    for (uint j = 0; j < M.n; j++) {
        for (uint k = 0; k < M.ku; k++) {
            if (j + k >= M.ku && j + k < M.n + M.ku) {
                if (std::abs(M.col[j][k]) > CInfo[5]) {
                    CInfo[5] = std::abs(M.col[j][k]);
                }
            }
        }
    }

    return CInfo;
}

void writeMatrixBandCtable(std::ofstream& output_file, const std::vector<double>& Infos) {
    output_file << "# c" << std::endl;
    output_file << std::setw(15) << "min(j<i) : " << std::fixed << std::setprecision(6) << Infos[0] << std::endl;
    output_file << std::setw(15) << "max(j<i) : " << std::fixed << std::setprecision(6) << Infos[1] << std::endl;
    output_file << std::setw(15) << "min(j=i) : " << std::fixed << std::setprecision(6) << Infos[2] << std::endl;
    output_file << std::setw(15) << "max(j=i) : " << std::fixed << std::setprecision(6) << Infos[3] << std::endl;
    output_file << std::setw(15) << "min(j>i) : " << std::fixed << std::setprecision(6) << Infos[4] << std::endl;
    output_file << std::setw(15) << "max(j>i) : " << std::fixed << std::setprecision(6) << Infos[5] << std::endl;
    output_file << std::endl;
}


//Write the Frobeinus ans the Infinite Normes
std::vector<double> ComputeDInfos(const mat_band& M) {
    std::vector<double> DInfos(2);

    DInfos[0] = normFrobeniusMatBand(M);
    DInfos[1] = normInftyMatBand(M);

    return DInfos;
}

void writeMatrixBandDtable(std::ofstream& output_file, const std::vector<double>& Infos) {
    output_file << "# d" << std::endl;
    output_file << std::setw(15) << "Frobenius: " << std::fixed << std::setprecision(6) << Infos[0] << std::endl;
    output_file << std::setw(15) << "Infinite: " << std::fixed << std::setprecision(6) << Infos[1] << std::endl;
    output_file << std::endl;
}


//Compute th LU and Cholesky factorization information
std::vector<double> ComputeEInfos(const mat_band& M) {
    std::vector<double> EInfos(10);
    clock_t c_start, c_end; // Variables used to measure the CPU time
    double t_LU, t_Cholesky; // LU and Cholesky decomposition time
    int check;

    mat_band M_LU = M;
    mat_band M_CHKY = M;

    const double nan = std::numeric_limits<double>::quiet_NaN(); // méthode non applicable
    vec X(M.n, 1.);
    vec B = matBandVec(M, X);

    // Measure the time taken for LU decomposition
    std::vector<uint> piv;
    c_start = clock(); int checkLU = LUBand(M_LU, piv); c_end = clock();
    t_LU = static_cast<double>(c_end - c_start) / CLOCKS_PER_SEC;

    // Solve the system using LU decomposition
    vec Xapp = checkLU == 0 ? LUBand_solve(M_LU, piv, B) : vec(M.n, nan);

    vec E(M.n);
    for (unsigned int i = 0u; i < M.n; ++i)
        E[i] = X[i] - Xapp[i];

    vec R = residualBand(B, M, Xapp);

    // Fill the EInfos vector with relevant information for LU decomposition
    EInfos[0] = t_LU;
    EInfos[1] = checkLU == 0 ? LUBand_det(M_LU, piv) : 0.;
    EInfos[2] = checkLU == 0 ? logAbsDiag(M_LU) : -std::numeric_limits<double>::infinity();
    EInfos[3] = normInftyVec(E);
    EInfos[4] = normInftyVec(R);

    // Measure the time taken for Cholesky decomposition
    c_start = clock(); check = CholeskiBand(M_CHKY); c_end = clock();
    t_Cholesky = static_cast<double>(c_end - c_start) / CLOCKS_PER_SEC;

    // Solve the system using Cholesky decomposition
    vec Xapp2 = check == 0 ? CholeskiBand_solve(M_CHKY, B) : vec(M.n, nan);
    for (unsigned int i = 0u; i < M.n; ++i)
        E[i] = X[i] - Xapp2[i];

    vec R2 = residualBand(B, M, Xapp2);

    //Fill the EInfos vector with relevant information for Cholesky decomposition
    if (check == -1) {
        EInfos[5] = EInfos[6] = EInfos[7] = nan;
    } else {
        EInfos[5] = t_Cholesky;
        EInfos[6] = CholeskiBand_det(M_CHKY);
        EInfos[7] = 2. * logAbsDiag(M_CHKY);
    }

    EInfos[3] = checkLU == 0 ? EInfos[3] : nan;
    EInfos[4] = checkLU == 0 ? EInfos[4] : nan;
    EInfos[8] = check == 0 ? normInftyVec(E) : nan;
    EInfos[9] = check == 0 ? normInftyVec(R2) : nan;

    return EInfos;
}

void writeMatrixBandEtable(std::ofstream& output_file, const std::vector<double>& Infos) {
    output_file <<  "# e" << std::endl;
    output_file << std::setw(35) << "LU factorization time: " << std::fixed << std::setprecision(6) << Infos[0] << std::endl;
    output_file << std::setw(35) << "LU determinant: " << std::scientific << std::setprecision(6) << Infos[1] << std::endl;
    output_file << std::setw(35) << "log(LU determinant): " << std::fixed << std::setprecision(6) << Infos[2] << std::endl;
    output_file << std::setw(35) << "LU error: " << std::scientific << std::setprecision(6) << Infos[3] << std::endl;
    output_file << std::setw(35) << "LU residual: " << std::scientific << std::setprecision(6) << Infos[4] << std::endl;

    output_file << std::setw(35) << "Cholesky factorization time: " << std::fixed << std::setprecision(6) << Infos[5] << std::endl;
    output_file << std::setw(35) << "Cholesky determinant: " << std::scientific << Infos[6] << std::endl;
    output_file << std::setw(35) << "log(Cholesky determinant): " << std::fixed << std::setprecision(6) << Infos[7] << std::endl;
    output_file << std::setw(35) << "Cholesky error: " << std::scientific << std::setprecision(6) << Infos[8] << std::endl;
    output_file << std::setw(35) << "Cholesky residual: " << std::scientific << std::setprecision(6) << Infos[9] << std::endl;

    output_file << std::endl;
}


//Compute the SOR factorization information
std::vector<double> ComputeFInfos(const mat_band& M) {
    std::vector<double> FInfos(6), Xapp(M.n, 0), X(M.n, 1);
    vec B = matBandVec(M, X);
    mat_band M_SOR = M;

    clock_t c_start, c_end; // Variables used to measure the CPU time
    double t_SOR; // SOR decomposition time

    double omega = 1, tolerance = 1e-6, max_iter = 10000, NbIteration;

    // Measure the time taken for SOR decomposition
    c_start = clock();
    NbIteration = SOR(M_SOR, B, Xapp, omega, max_iter, tolerance);
    c_end = clock();

    vec E(M.n);
    for (unsigned int i = 0u; i < M.n; ++i)
        E[i] = X[i] - Xapp[i];

    vec R = residualBand(B, M, Xapp);

    // If the number of iterations is equal to the maximum number of iterations, the SOR decomposition failed
    if (NbIteration == max_iter) {
        t_SOR = 0;
    } else {
        t_SOR = static_cast<double>(c_end - c_start) / CLOCKS_PER_SEC;
    }

    FInfos[0] = omega;
    FInfos[1] = tolerance;
    FInfos[2] = t_SOR;
    FInfos[3] = NbIteration;
    FInfos[4] = normInftyVec(E);
    FInfos[5] = normInftyVec(R);

    return FInfos;
}


void writeMatrixBandFtable(std::ofstream& output_file, const std::vector<double>& Infos) {
    output_file << "# f" << std::endl;
    output_file << std::setw(15) << "omega: " << std::fixed << std::setprecision(6) << Infos[0] << std::endl;
    output_file << std::setw(15) << "tolerance: " << std::scientific << std::setprecision(6) << Infos[1] << std::endl;
    output_file << std::setw(35) << "SOR factorization time: " << std::fixed << std::setprecision(6) << Infos[2] << std::endl;
    output_file << std::setw(35) << "Number of iterations: " << Infos[3] << std::endl;
    output_file << std::setw(35) << "SOR error: " << std::scientific << std::setprecision(6) << Infos[4] << std::endl;
    output_file << std::setw(35) << "SOR residual: " << std::scientific << std::setprecision(6) << Infos[5] << std::endl;

    output_file << std::endl;
}