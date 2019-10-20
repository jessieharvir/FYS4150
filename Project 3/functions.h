#ifndef FUNC_H
#define FUNC_H
#include <iostream>
#include <armadillo>
#include <cmath>
#include <string>
#include "time.h"
#include <fstream>

using namespace std;
using namespace arma;

class functions{
public:
    functions();
    double func_cartesian(double x1, double y1, double z1, double x2, double y2, double z2);
    double func_polar(double r1, double theta1, double phi1, double r2, double theta2, double phi2);
    double gammln( double xx);
    void gauss_laguerre(double *x, double *w, int n, double alf);
    void gauss_legendre(double x1, double x2, double x[], double w[], int n);
    void loopyloop_cartesian(int n, double a, double b);
    void loopyloop_polar(int n);
    void monte_carlo_cartesian(int n, double a, double b);
    void monte_carlo_polar(int n);
    void MPI_monte_carlo_cartesian(int n, double a, double b);
    void MPI_monte_carlo_polar(int n);
    void tester_bestfunc();
    void tester_computed_exact_values();
    };

#endif //FUNC_H