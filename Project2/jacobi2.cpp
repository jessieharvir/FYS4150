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

int jacobi::jacobii(int n, int interact, double conv, double wr,mat& a, vec& r, mat& v) {
    cout.precision(5);
    double aip=0, aiq=0, vpi=0, vqi=0;
    double tau=0, t=0, s=0, c=0;//tan(theta), sin(theta), cos(theta)    
    int count=1;                //count of iterations
    int count_old=count-10;     //keep track of every 10th iteration    
    int p=n-1, q=n-2;           //off diag all same value to start
                                //pick last as first maximum
    clock_t start, end;

    if(n<=10){
        cout<<"Before diagonalization"<<endl;
        print_vals(a,v,n,conv);
        cout<<endl;
    }

    double app=a(p,p);
    double aqq=a(q,q);
    double apq=a(p,q);
    
    start=clock();
    
    while(abs(apq)>conv){
        if(count>1){
            apq=0;
            find_max(a,p,q,apq,n);
        }

        //calculate sin(theta) and cos(theta)
        aqq=a(q,q);
        app=a(p,p);
        tau=(aqq-app)/(2*apq);
        if(tau>0)
            t=1/(tau+sqrt(1+tau*tau));
        else
            t=-1/(-tau+sqrt(1+tau*tau));   
        c=1/sqrt(1+t*t);
        s=c*t;

        //calculate new matrix elements and vectors
        for(int i=0;i<n;i++){
            if(i!=p && i!=q){
                aip=a(i,p);
                aiq=a(i,q);
                a(i,p)=aip*c-aiq*s;
                a(p,i)=aip*c-aiq*s;
                a(i,q)=aiq*c+aip*s;
                a(q,i)=aiq*c+aip*s;
            }
            //vpi=v(p,i);
            //vqi=v(q,i);
            vpi=v(i,p);
            vqi=v(i,q);
            //v(p,i)=c*vpi-s*vqi;
           // v(q,i)=c*vqi+s*vpi;
            v(i,p)=c*vpi-s*vqi;
            v(i,q)=c*vqi+s*vpi;
        }
        a(p,p)=app*c*c-2*apq*c*s+aqq*s*s;
        a(q,q)=app*s*s+2*apq*c*s+aqq*c*c;
        a(p,q)=0;
        a(q,p)=0;
        
        count++;
    }
    
    end=clock();
    
    if(n<=10){
        cout<<"After diagonalization"<<endl;
        print_vals(a,v,n,conv);
        cout<<endl;
    }

    cout<<"Diagonalization took "<<count<<" iterations"<<endl;
    cout<<scientific<<"CPU time (sec) : "<<((double)end-(double)start)/CLOCKS_PER_SEC<<endl;
    
    return 0;
}


//get first three eigenvectors
mat jacobi::get_eigenvecs(mat a, mat v, int n){
    vector<double>eigenvals=get_eigenvals(a,n);
    mat vecs(3,n);
    for(int i=0;i<3;i++){
        for(int j=0;j<n;j++){
            if(a(j,j)==eigenvals[i]){
                for(int k=0;k<n;k++){
                      vecs(i,k)=v(k,j);
                }
             }
         }
    }
    return vecs;
}

//get eigenvalues in order
vector<double> jacobi::get_eigenvals(mat a,int n){
    vector<double>eigen;
    for(int i=0;i<n;i++){
        eigen.push_back(a(i,i));
    }
    sort (eigen.begin(), eigen.begin()+n);
    return eigen;
}

//initialize matrix/vectors
void jacobi::initialize(int n, double h, mat& a, vec& r, mat& v,int interact,double wr){
    //initialize x values
    r(0)=h;
    for (int i=1; i<n ;i++){
        r(i)=r(i-1)+h;
    }
    
    //initialize matrix and vector
    for (int i=0;i<n;i++){
        for (int j=0;j<n;j++){
            if(i==j && interact==0){
                a(i,j)=2/(h*h)+r(i)*r(i);
                v(i,j)=1;
            }
            else if (i==j && interact==1){
                a(i,j)=2/(h*h)+wr*wr*r(i)*r(i)+1/r(i);
                v(i,j)=1;
            }
            else if (i==j+1 or i==j-1){
                a(i,j)=-1/(h*h);
            } 
            else{
                a(i,j)=0;
                v(i,j)=0;
            }
        }
    }
}

//find maximum non-diag matrix elements
void jacobi::find_max(mat a,int& p,int& q,double& apq,int n){
    for (int i=0;i<n;i++){
         for (int j=0;j<n;j++){
            if(i!=j && abs(a(i,j))>=abs(apq)){
                apq=a(i,j);
                p=i;
                q=j;
            }
         }
    }
}

//print matrix and eigenvectors
void jacobi::print_vals(mat A, mat v,int n,double conv){
    cout<<"A: ";
    for (int i=0;i<n;i++){
        if(i>0){
            cout<<"   ";
         }
        for (int j=0;j<n;j++){
            if(abs(A(i,j))>conv)
                cout<<fixed<<A(i,j)<<" ";
            else cout<<"0.000 ";
        }
        cout<<endl;
    }
    for (int i=0;i<n;i++){
        cout<<"v"<<i<<": ";
        for (int j=0;j<n;j++){
            if(abs(v(j,i))>conv)
                cout<<fixed<<v(j,i)<<" ";
            else cout<<"0.000 ";
        }
        cout<<endl;
    }  
}