#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "jacobi.h"
#include "armadillo"
jacobi jacobi_method;

TEST_CASE("Testing creation of diagonal matrix and max A(i,j)"){
    int n = 4;
    mat A(n,n,fill::zeros);

    double d = 2.0; 
    double a = -1.0;
    int p, q; 
    A = jacobi_method.create_matrix(A, n, d, a, 0.0, 0.0, 0.0);
    REQUIRE(A(2,2)==2);
    REQUIRE(A(1,2)==-1);
    REQUIRE(A(3,0)==0);

    A(1,3) = 5.6;
    double maxnondiag = jacobi_method.offdiag(A, p, q, n);
    REQUIRE(maxnondiag==Approx(5.6));
    
}

TEST_CASE("Testing eigenvalues calculated by Jacobis method"){
    int n = 4;
    double tolerance = 1.0E-15;
    int iterations = 0;
    int maxiter = 100000;
    double rho0 = 0;
    double rho_max = 4.5;
    double d = 2.0; 
    double a = -1.0;
    mat A(n,n,fill::zeros);
    mat R(n,n,fill::eye);
    A = jacobi_method.create_matrix(A, n, d, a, rho_max, false, 0.0);
    double maxnondiag = tolerance*2;
    vec eig_arma = sort(eig_sym(A));
    while (maxnondiag > tolerance && iterations <= maxiter){
        int p, q;
        maxnondiag = jacobi_method.offdiag(A, p, q, n);
        jacobi_method.jacobi_rotate(A, R, p, q, n);
        iterations++;
    }
    vec eig_jacobi = sort(diagvec(A));
    REQUIRE(eig_arma[0]==Approx(eig_jacobi[0]));

}
/*


TEST_CASE("Testing eigenvalues"){
    int n=4,interact=0;
    double conv=0.001,wr=0.01, pmin=0, pmax=10,h = (pmax-pmin)/(double(n));
    n=n-1;
    mat a = zeros<mat>(n,n);
    mat v = zeros<mat>(n,n);
    vec r(n);
    //initialize matrices and vector
    jacobi_method.initialize(n,h,a,r,v,0,0);
    //do jacobi algorithm until convergence
    jacobi_method.jacobii(n,interact,conv,wr,a,r,v);
    //get eigenvalue vector
    vector<double>eigen=jacobi_method.get_eigenvals(a,n);
    
    REQUIRE(eigen[0]==Approx(6.56863));
    REQUIRE(eigen[1]==Approx(25.32055));
    REQUIRE(eigen[2]==Approx(56.57082));
}
TEST_CASE("Testing eigenvector orthogonality"){
    int n=4,interact=0;
    double conv=0.001,wr=0.01, pmin=0, pmax=10,h = (pmax-pmin)/(double(n));
    n=n-1;
    mat a = zeros<mat>(n,n);
    mat v = zeros<mat>(n,n);
    vec r(n);
    //initialize matrices and vector
    jacobi_method.initialize(n,h,a,r,v,0,0);
    //do jacobi algorithm until convergence
    jacobi_method.jacobii(n,interact,conv,wr,a,r,v);
    mat eigenvec=jacobi_method.get_eigenvecs(a,v,n);

    //test eigen vectors
    REQUIRE(eigenvec(0,0)==Approx(0.9996).epsilon(0.001));
    REQUIRE(eigenvec(0,1)==Approx(0.00853));
    REQUIRE(eigenvec(0,2)==Approx(0.0000).epsilon(0.001));
    REQUIRE(eigenvec(1,0)==Approx(-0.00853));
    REQUIRE(eigenvec(1,1)==Approx(0.99995).epsilon(0.001));
    REQUIRE(eigenvec(1,2)==Approx(0.00512).epsilon(0.001));
    REQUIRE(eigenvec(2,0)==Approx(0.0000).epsilon(0.001));
    REQUIRE(eigenvec(2,1)==Approx(-0.00512).epsilon(0.001));
    REQUIRE(eigenvec(2,2)==Approx(0.9999).epsilon(0.001));
    
    //test eigen vector orthogonality
    //dot1=v0*v1=0
    double dot1=eigenvec(0,0)*eigenvec(1,0)+eigenvec(0,1)*eigenvec(1,1)
        +eigenvec(0,2)*eigenvec(1,2);
    //dot2=v0*v0=1
    double dot2=eigenvec(0,0)*eigenvec(0,0)+eigenvec(0,1)*eigenvec(0,1)
        +eigenvec(0,2)*eigenvec(0,2);
    REQUIRE(dot1==Approx(0.000).epsilon(0.01));
    REQUIRE(dot2==Approx(1.000));
}
*/