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
    int jacobii(int n, int interact, double conv, double wr,mat& a, vec& r, mat& v);
    mat get_eigenvecs(mat a, mat v, int n);
    vector<double> get_eigenvals(mat a,int n);
    void initialize(int n, double h, mat& a, vec& r, mat& v,int interact,double wr);
    void find_max(mat a,int& p,int& q,double& apq,int n);
    void print_vals(mat A, mat v,int n,double conv);
};

#endif //JACOBI_H