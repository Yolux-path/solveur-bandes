#include <iostream>  // The std library for i/o
#include <fstream>   // ... for file streams
#include <sstream>   // ... for string streams
#include <iomanip>   // Manipulate i/o, e.g. formatting

#include <algorithm>
#include <cmath>     // The functions power and square root
#include <limits>    // Limits in FP numbers (IEEE-754)

#include <vector>    // The std library for the vector container type
#include "mat_vec_norm.h" 

using namespace std; // to avoid writing std:: everywhere, this is not good
                     // practise

typedef unsigned int uint;  // Use uint for unsigned int


// Computes 1-norm of vector
double norm1Vec(const vec& B){
    uint n=B.size();
    double Norm1 =0;
    for(uint i=0; i<n;i++){
        Norm1 += fabs(B[i]);
    }
    return Norm1;
}

// Computes 2-norm of vector
double norm2Vec(const vec& B) {
    uint n = B.size();
    double Norm2=0;
    for(uint i=0;i<n;i++){
        Norm2+=pow(B[i],2);
    }
    return sqrt(Norm2);
}

double normFrobeniusMat(const matrix& M){
    double NormFrobe=0;
    
    for(uint i=0;i<M.n;i++){
        for(uint j=0;j<M.m;j++){
            NormFrobe+=pow(M.a[i][j],2);
        }
    }

    return sqrt(NormFrobe);
}

double normFrobeniusMatBand(const mat_band& M){
    
    double NormFrobe=0;
    
    for(uint j = 0; j < M.n; j++){
        for(uint k = 0; k < M.ku+M.kl+1 ; k++){
            if(j+k >= M.ku && j+k < M.n+M.ku){
                
                NormFrobe += pow( abs(M.col[j][k]) , 2 );
            }
        }
    }

    return sqrt(NormFrobe);
}

// Computes infinite norm of vector
double normInftyVec(const vec& B) {
 
    double res = 0.e0;
    for (uint i=0u; i<B.size(); i++) {
        if (std::isnan(B[i])) return B[i]; // sinon un NaN passerait pour une erreur nulle
        if (fabs(B[i])>res) res = fabs(B[i]);
    }
    return res;
}

// Computes infinite norm of matrix (max over rows of sum of abs of coefficients)
double normInftyMat(const matrix& M){
    double NormInf=0;
    vec SumParL(M.n,0.);
    for(uint i=0;i<M.n;i++){
        for(uint j=0;j<M.m;j++){
            SumParL[i]+=abs(M.a[i][j]);
        }
    }
    NormInf =*max_element(SumParL.begin(),SumParL.end());
    return NormInf;
}

double normInftyMatBand(const mat_band& M){

    double NormInf = 0;
    vec SumparL(M.n, 0); 

    for(uint j = 0; j < M.n; j++){
        for(uint k = 0; k < M.ku+M.kl+1 ; k++){

            if(j+k >= M.ku && j+k < M.n+M.ku){
                SumparL[k + j - M.ku ] += abs(M.col[j][k]) ;
            }

        }
    }

    NormInf =*max_element(SumparL.begin(),SumparL.end());
    return NormInf;
            

}