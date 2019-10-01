#include <iostream>
#include <armadillo>
#include <cmath>
#include <string>
#include "time.h"
#include <fstream> 
#include "jacobi.h"
#include <chrono>

using namespace std;
using namespace arma;

mat create_matrix(mat M, int n, double d, double a);
void test_diagonals(mat A, mat R, int n, double tolerance, int iterations, int maxiter, double a, double d, double rho_max);
void quantumdots_3D_oneelec(mat A, mat R, int n, double tolerance, int iterations, int maxiter, double a, double d, double rho_max);
void quantumdots_3D_twoelec(mat A, mat R, int n, double tolerance, int iterations, int maxiter, double a, double d, double rho_max, double omega);

int main(){
    int n = 200;
    mat A(n,n,fill::zeros);
    mat R(n,n,fill::eye);
    double tolerance = 1.0E-15;
    int iterations = 0;
    int maxiter = 100000;
    double rho0 = 0;
    double rho_max = 4.5;
    //double h = (rho_max-rho0)/(n+1);
    double d = 2;
    double a = -1;

    //test_diagonals(A, R, n, tolerance, iterations, maxiter, a, d, rho_max);
    /*
    vec rhos = linspace<vec>(4.4, 4.5, 10); 
    for (int i=1; i<10; i++){
        cout << "rho = " << rhos[i] << endl;
        quantumdots_3D_oneelec(A, R, n, tolerance, iterations, maxiter, a, d, rhos[i]);
    } */  
    
    test_diagonals(A, R, 5, tolerance, iterations, maxiter, a, d, rho_max);
    test_diagonals(A, R, 20, tolerance, iterations, maxiter, a, d, rho_max);
    test_diagonals(A, R, 100, tolerance, iterations, maxiter, a, d, rho_max);
    test_diagonals(A, R, 200, tolerance, iterations, maxiter, a, d, rho_max);
    
    /*
    quantumdots_3D_twoelec(A, R, n, tolerance, iterations, maxiter, a, d, rho_max, 0.01);
    quantumdots_3D_twoelec(A, R, n, tolerance, iterations, maxiter, a, d, rho_max, 0.5);
    quantumdots_3D_twoelec(A, R, n, tolerance, iterations, maxiter, a, d, rho_max, 1.0);
    quantumdots_3D_twoelec(A, R, n, tolerance, iterations, maxiter, a, d, rho_max, 5.0);
    */
}



void test_diagonals(mat A, mat R, int n, double tolerance, int iterations, int maxiter, double a, double d, double rho_max){
    jacobi jacobi_method;
    A = jacobi_method.create_matrix(A, n, d, a, rho_max, false, 0.0);
    double maxnondiag = tolerance*2;
    auto t1a = std::chrono::high_resolution_clock::now();
    vec eig = eig_sym(A);
    //cout << sort(eig) << endl;
    auto t2a = std::chrono::high_resolution_clock::now();
    auto durationa = std::chrono::duration_cast<std::chrono::microseconds>( t2a - t1a ).count();
    cout << "Time used for n = " << n << " to find eigenvalues with armadillo " << durationa*1e-6 << " sec." << endl; 
    
      
    auto t1 = std::chrono::high_resolution_clock::now();
    while (maxnondiag > tolerance && iterations <= maxiter){
        int p, q;
        maxnondiag = jacobi_method.offdiag(A, p, q, n);
        jacobi_method.jacobi_rotate(A, R, p, q, n);
        iterations++;
    }
    cout << "Number of iterations: " << iterations << endl;
    //vec new_eigenvec = diagvec(A);
    //cout << sort(new_eigenvec) << endl;
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    cout << "Time used for n = " << n << " to do the jacobi rotate " << duration*1e-6 << " sec." << endl;   
}

void quantumdots_3D_oneelec(mat A, mat R, int n, double tolerance, int iterations, int maxiter, double a, double d, double rho_max){
    jacobi jacobi_method;
    A = jacobi_method.create_matrix(A, n, d, a, rho_max, true, 0.0);
    double maxnondiag = tolerance*2;
    auto t1 = std::chrono::high_resolution_clock::now();
    while (maxnondiag > tolerance && iterations <= maxiter){
        int p, q;
        maxnondiag = jacobi_method.offdiag(A, p, q, n);
        jacobi_method.jacobi_rotate(A, R, p, q, n);
        iterations++;
    }  
    cout << "Number of iterations: " << iterations << endl;
    vec eig_sorted = sort(diagvec(A));
    cout << "Eigenvalues: " << eig_sorted[0] << " " << eig_sorted[1] << " " << eig_sorted[2] << " " << eig_sorted[3] << " " << endl;
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    cout << "Time used for n = " << n << " to create matrix A and solving eigenvalues for one electron " << duration*1e-6 << " sec." << endl;  
}

void quantumdots_3D_twoelec(mat A, mat R, int n, double tolerance, int iterations, int maxiter, double a, double d, double rho_max, double omega){
    jacobi jacobi_method;

    A = jacobi_method.create_matrix(A, n, d, a, rho_max, true, omega);
    double maxnondiag = tolerance*2;
    auto t1 = std::chrono::high_resolution_clock::now();
    while (maxnondiag > tolerance && iterations <= maxiter){
        int p, q;
        maxnondiag = jacobi_method.offdiag(A, p, q, n);
        jacobi_method.jacobi_rotate(A, R, p, q, n);
        iterations++;
    }
    cout << "Omega_r = " << omega << endl;
    cout << "Number of iterations: " << iterations << endl;
    vec eig_sorted = sort(diagvec(A));
    cout << "Eigenvalues: " << eig_sorted[0] << " " << eig_sorted[1] << " " << eig_sorted[2] << " " << eig_sorted[3] << " " << endl;  
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    cout << "Time used for n = " << n << " to create matrix A and solving eigenvalues for one electron " << duration*1e-6 << " sec." << endl;  
}