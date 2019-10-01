#ifndef JACOBI_H
#define JACOBI_H
#include <iostream>
#include <armadillo>
#include <cmath>
#include <string>
#include "time.h"
#include <fstream>

using namespace std;
using namespace arma;

class jacobi{
public:
    jacobi();
    double offdiag(mat& A, int& p, int& q, int n);
    void jacobi_rotate(mat &A, mat &R, int k, int l, int n);
    mat create_matrix(mat M, int n, double d, double a, double rho_max, bool potential, double omega);
};

#endif //JACOBI_H