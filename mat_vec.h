#ifndef __MAT_VEC_H_INCLUDED__
#define __MAT_VEC_H_INCLUDED__ 

#include <vector> 

typedef std::vector<double> vec;

// A struct to store a matrix
struct matrix {
  unsigned int n,m; 
  std::vector<std::vector<double>> a;
};

// A struct to store a band matrix
struct mat_band {
  unsigned int n;
  unsigned int kl,ku;
  std::vector<std::vector<double>> col;
  
};

#endif 
