#include "jacobi.h"

jacobi::jacobi(){}
//  the offdiag function, using Armadillo
double jacobi::offdiag(mat& A, int& p, int& q, int n){
   double maxnondiag = 0;
   for (int i = 0; i < n; ++i)
   {
       for ( int j = i+1; j < n; ++j)
       {
            double aij = fabs(A(i,j));
           if ( aij > maxnondiag)
           { 
              maxnondiag = aij;  p = i; q = j;
           }
       }
   }
   return maxnondiag;
}
// more statements

void jacobi::jacobi_rotate (mat &A, mat &R, int k, int l, int n ){
  double s, c;
  if ( A(k,l) != 0.0 ) {
    double t, tau;
    tau = (A(l,l) - A(k,k))/(2 * A(k,l));
    
    if ( tau >= 0 ) {
      t = 1.0/(tau + sqrt(1.0 + tau*tau));
    }
    else {
      t = -1.0/(-tau +sqrt(1.0 + tau*tau));
    }
    c = 1/sqrt(1+t*t);
    s = c*t;
    } 
    else {
    c = 1.0;
    s = 0.0;
    }
    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    a_kk = A(k,k);
    a_ll = A(l,l);
    A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
    A(l,l) = s*s*a_kk + 2.0*c*s*A(k,l) + c*c*a_ll;
    A(k,l) = 0.0;  // hard-coding non-diagonal elements by hand
    A(l,k) = 0.0;  // same here
    for ( int i = 0; i < n; i++ ) {
        if ( i != k && i != l ) {
            a_ik = A(i,k);
            a_il = A(i,l);
            A(i,k) = c*a_ik - s*a_il;
            A(k,i) = A(i,k);
            A(i,l) = c*a_il + s*a_ik;
            A(l,i) = A(i,l);
    }
//  And finally the new eigenvectors
    r_ik = R(i,k);
    r_il = R(i,l);

    R(i,k) = c*r_ik - s*r_il;
    R(i,l) = c*r_il + s*r_ik;
  }
  return;
} // end of function jacobi_rotate

mat jacobi::create_matrix(mat M, int n, double d, double a, double rho_max, bool potential, double omega){
  vec rho = linspace<vec>(0, rho_max, n+2);
  double h = rho_max/(n+1);
  a = a/(h*h);
  d = d/(h*h);
  M(0,0) = d;
  M(n-1,n-1) = d;
  M(0,1) = a;
  M(1,0) = a;
  M(n-1,n-2) = a;
  M(n-2,n-1) = a;

  if (potential){
    if (omega != 0.0){
      M(0,0) = d + rho[1]*rho[1] * omega*omega + 1/rho[1];
      M(n-1,n-1) = d + rho[n]*rho[n] * omega*omega + 1/rho[n];
      for (int i = 1; i <= n-2; i++){
          M(i,i) = d + rho[i+1]*rho[i+1] * omega*omega + 1/rho[i+1];
          M(i, i+1) = a;
          M(i, i-1) = a;
      }
    }
    else{ 
      M(0,0) = d + rho[1]*rho[1];
      M(n-1,n-1) = d + (n-1)*h*(n-1)*h;
      for (int i = 1; i <= n-2; i++){
        M(i,i) = d + rho[i+1]*rho[i+1];
        M(i, i+1) = a;
        M(i, i-1) = a;
      }
    }
  }
  else {
      for (int i = 1; i <= n-2; i++){
        M(i,i) = d;
        M(i, i+1) = a;
        M(i, i-1) = a;
      }
  }
  return M;
}