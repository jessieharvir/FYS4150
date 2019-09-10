#include <iostream>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <string>
#include "armadillo"

using namespace std;
using namespace arma;

ofstream ofile;

void prob1_b(int n);
void prob1_c(int n);
double exact_solution(double x);
double analytic_solution(double x);
void LU_decomp(int n);


int main()
{
    int n = 10;
    //cout << "Please provide with a n-value, either 10, 100 or 1000" << endl;
    //cin >> n;
    //prob1_b(n);

    int n_values[7] = {10, 100, 1000, 10000, 100000, 1000000, 10000000};
    for (int i = 0; i <= 7; i++){
        prob1_b(n_values[i]);
        prob1_c(n_values[i]);
        //LU_decomp(n_values[i]);
    }

    //LU_decomp(n);
    return 0;
}

double analytic_solution(double x){
    return 1.0-(1-exp(-10))*x-exp(-10*x);
}

double exact_solution(double x){
    return 100*exp(-10*x);
}

void prob1_b(int n){
    // Diagonals of matrix
    double *a = new double[n+2];
    double *b = new double[n+2];
    double *b_tilde = new double[n+2];
    double *c = new double[n+2];

    // Analytical and numerucal solution
    double *u = new double[n+2];
    double *v = new double[n+2];

    // Step size
    double h = 1.0/(n+1.0);

    // Steps
    double *x = new double[n+2];

    // Right hand side of equation
    double *f = new double[n+2];
    double *f_tilde = new double[n+2];

    // Relative error
    double *error = new double[n+2];

    auto t1 = std::chrono::high_resolution_clock::now();
    // the vectors
    for (int i = 0; i < n+1; i++){
        x[i] = i*h;

        a[i] = 2;
        b[i] = -1;
        c[i] = -1;

        f[i] = h*h*exact_solution(x[i]);
        u[i] = analytic_solution(x[i]);
    }
    // Gaussian elimination
    //forward subst
    f_tilde[0] = f[0];
    b_tilde[0] = b[0];
    for(int i=1; i<n+1; i++){
        b_tilde[i] = b[i] - a[i-1]*c[i-1]/b_tilde[i-1];
        f_tilde[i] = f[i] - a[i-1]*f_tilde[i-1]/b_tilde[i-1];
    }
    
    //backwards subst
    v[0] = v[n] = 0;
    v[n-1] = f_tilde[n-1]/b_tilde[n-1];
    for(int i=n-2; i>=1; i--){
        v[i] = (f_tilde[i]-c[i]*v[i+1])/b_tilde[i];
        error[i] = log10(abs((v[i]-u[i])/u[i]));

    }
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    cout << "Time used for n = " << n << " with the general gaussian elimination is " << duration*1e-6 << " sec." << endl;

    // write results to outputfile
    string name = "/home/jessie/Dokumenter/skole/H19/FYS4150/Project_1/Project1/prob1_b_n" + to_string(n) + ".txt";
    ofile.open(name);
    if (ofile.is_open()) {
        ofile << setiosflags(ios::showpoint | ios::uppercase);
        for (int i = 1; i < n;i++) {
             ofile << setw(15) << setprecision(8) << x[i];
             ofile << setw(15) << setprecision(8) << u[i];
             ofile << setw(15) << setprecision(8) << v[i];
             ofile << setw(15) << setprecision(8) << error[i] << endl;
        }
    }
    else cout << "Unable to open file";
    ofile.close();


    // delete arrays
    delete [] a; delete [] b; delete [] b_tilde; delete [] c; delete [] f; delete [] f_tilde; delete [] x; delete [] u; delete [] v;
}

void prob1_c(int n){
    double *b_tilde = new double[n+2];

    // Analytical and numerucal solution
    double *u = new double[n+2];
    double *v = new double[n+2];

    // Step size
    double h = 1.0/(n+1.0);

    // Steps
    double *x = new double[n+2];

    // Right hand side of equation
    double *f = new double[n+2];
    double *f_tilde = new double[n+2];

    // Relative error
    double *error = new double[n+2];

    auto t1 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < n+1; i++){
        x[i] = i*h;
        f[i] = h*h*exact_solution(x[i]);
        u[i] = analytic_solution(x[i]);
    }

    // Gaussian elimination
    //forward subst
    f_tilde[0] = f[0];
    b_tilde[0] = 2;
    for(int i=1; i<n+1; i++){
        b_tilde[i] = 2 - 1/b_tilde[i-1];
        f_tilde[i] = f[i] + f_tilde[i-1]/b_tilde[i-1];
    }

    //backwards subst
    v[0] = v[n] = 0;
    v[n-1] = f_tilde[n-1]/b_tilde[n-1];
    for(int i=n-2; i>=1; i--){
        v[i] = (f_tilde[i]+v[i+1])/b_tilde[i];
        error[i] = log10(abs((v[i]-u[i])/u[i]));
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    cout << "Time used for n = " << n << " with the special gaussian elimination is " << duration*1e-6 << " sec." << endl;

    // write results to outputfile
    string name = "/home/jessie/Dokumenter/skole/H19/FYS4150/Project_1/Project1/prob1_c_n" + to_string(n) + ".txt";
    ofile.open(name);
    if (ofile.is_open()) {
        ofile << setiosflags(ios::showpoint | ios::uppercase);
        for (int i = 1; i < n;i++) {
             ofile << setw(15) << setprecision(8) << x[i];
             ofile << setw(15) << setprecision(8) << u[i];
             ofile << setw(15) << setprecision(8) << v[i];
             ofile << setw(15) << setprecision(8) << error[i] << endl;
        }
    }
    else cout << "Unable to open file";
    ofile.close();


    // delete arrays
    delete [] b_tilde; delete [] f; delete [] f_tilde ; delete [] x; delete [] u; delete [] v;
}

void LU_decomp(int n){
    // Exact solution
    vec v = zeros(n);


    // Right hand side of equation
    vec f = zeros<vec>(n);

    // Step size
    double h = 1.0/(n+1.0);

    // Steps
    double *x = new double[n];

    auto t1 = std::chrono::high_resolution_clock::now();
    // the vectors
    for (int i = 0; i < n-1; i++){
        x[i] = i*h;
        f(i) = h*h*exact_solution(x[i]);
    }


    mat A = zeros<mat>(n,n);
    //cout << A << endl;

    A(0,0) = 2.;
    A(n-1,n-1) = 2;
    A(0,1) = -1;
    A(1,0) = -1;
    A(n-1,n-2) = -1;
    A(n-2,n-1) = -1;
    for (int i = 1; i <= n-2; i++){
        A(i,i) = 2;
        A(i, i+1) = -1;
        A(i, i-1) = -1;
    }
    mat L, U;
    lu(L, U, A);
    vec y = solve(trimatl(L), f);
    v = solve(trimatu(U), y);
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    cout << "Time used for n = " << n << " with the LU-decomposition is " << duration*1e-6 << " sec." << endl;

}
