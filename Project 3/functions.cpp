#include <cmath>
#include <mpi.h>
#include "functions.h"
#include <cassert>


functions::functions(){}
double functions::func_cartesian(double x1, double y1, double z1, double x2, double y2, double z2){
	if  ((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2) > 1.0E-8)
		return exp(-4*(sqrt(x1*x1+y1*y1+z1*z1)+sqrt(x2*x2+y2*y2+z2*z2))) 
		          / sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
	else 
		return 0;
}
double functions::func_polar(double r1, double theta1, double phi1, double r2, double theta2, double phi2){
	double cosb = cos(theta1)*cos(theta2) + sin(theta1)*sin(theta2)*cos(phi1-phi2);
	double f = exp(-3*(r1+r2))*r1*r1*r2*r2*sin(theta1)*sin(theta2)/sqrt(r1*r1+r2*r2-2*r1*r2*cosb);
	if(r1*r1+r2*r2-2*r1*r2*cosb > 1.0E-8)
		return f;
	else 
		return 0;
}
double functions::gammln( double xx){
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}
void functions::gauss_laguerre(double *x, double *w, int n, double alf){
	int i,its,j;
	double ai;
	double p1,p2,p3,pp,z,z1;
    int MAXIT = 10;
    double EPS = 3.0e-14;

	for (i=1;i<=n;i++) {
		if (i == 1) {
			z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
		} else if (i == 2) {
			z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
		} else {
			ai=i-2;
			z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
				(1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
		}
		for (its=1;its<=MAXIT;its++) {
			p1=1.0;
			p2=0.0;
			for (j=1;j<=n;j++) {
				p3=p2;
				p2=p1;
				p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
			}
			pp=(n*p1-(n+alf)*p2)/z;
			z1=z;
			z=z1-p1/pp;
			if (fabs(z-z1) <= EPS) break;
		}
		if (its > MAXIT) cout << "too many iterations in gaulag" << endl;
		x[i]=z;
		w[i] = -exp(gammln(alf+n)-gammln((double)n))/(pp*n*p2);
	}
}
void functions::gauss_legendre(double x1, double x2, double x[], double w[], int n){
    int         m,j,i;
    double      z1,z,xm,xl,pp,p3,p2,p1;
    double      const  pi = 3.14159265359; 
    double      *x_low, *x_high, *w_low, *w_high;
    #define   ZERO       1.0E-10

    m  = (n + 1)/2;                             // roots are symmetric in the interval
    xm = 0.5 * (x2 + x1);
    xl = 0.5 * (x2 - x1);

    x_low  = x;                                       // pointer initialization
    x_high = x + n - 1;
    w_low  = w;
    w_high = w + n - 1;

    for(i = 1; i <= m; i++) {                             // loops over desired roots
        z = cos(pi * (i - 0.25)/(n + 0.5));
        do {
            p1 =1.0;
	        p2 =0.0;

	        for(j = 1; j <= n; j++) {
	            p3 = p2;
	            p2 = p1;
	            p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3)/j;
            }
	        pp = n * (z * p1 - p2)/(z * z - 1.0);
	        z1 = z;
	        z  = z1 - p1/pp;                   // Newton's method
        } 
        while(fabs(z - z1) > ZERO);
            *(x_low++)  = xm - xl * z;
            *(x_high--) = xm + xl * z;
            *w_low      = 2.0 * xl/((1.0 - z * z) * pp * pp);
            *(w_high--) = *(w_low++);
    }
}
void functions::loopyloop_cartesian(int n, double a, double b){
	double *roots = new double[n];
	double *weights = new double[n];
    auto t1 = std::chrono::high_resolution_clock::now();
	gauss_legendre(a, b, roots, weights, n);

	double integral = 0;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < n; k++) {
				for (int l = 0; l < n; l++) {
					for (int m = 0; m < n; m++) {
						for (int o = 0; o < n; o++) {
							integral += weights[i]*weights[j]*weights[k]*weights[l]*weights[m]*weights[o] * func_cartesian(roots[i],roots[j],roots[k],roots[l],roots[m],roots[o]);
						}
					}
				}
			}
		}
	}
    double const pi = 3.14159265359;
    double exact = 5*pi*pi/(16.0*16.0);
    cout << "N          = " << n << endl;
    cout << "Numeric    = " << integral << endl;
    cout << "Exact      = " << exact << endl;
    cout << "Rel. error = " << (abs(exact-integral))/exact << endl;
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    cout << "Time used for brute force Cartesian: " << duration*1e-6 << " sec." << endl; 
}
void functions::loopyloop_polar(int n){
	double const pi = 3.14159265359;
	double *theta = new double[n];
	double *phi = new double[n];
	double*r = new double[n];
	
	double *weight_theta = new double[n];
	double *weight_phi = new double[n];
	double *weight_r = new double[n];
    auto t1 = std::chrono::high_resolution_clock::now();
	gauss_legendre(0, pi, theta, weight_theta, n);
	gauss_legendre(0, 2*pi, phi, weight_phi, n);

	gauss_laguerre(r, weight_r, n, 0);

	double theta1, theta2, dtheta1, dtheta2;
	double phi1, phi2, dphi1, dphi2;
	double r1, r2, dr1, dr2;
	double integral = 0;

	for (int i = 0; i < n; i++) {
		dr1 = weight_r[i+1];
		r1 = r[i+1];
		for (int j = 0; j < n; j++) {
			dr2 = weight_r[j+1];
			r2 = r[j+1];
			for (int k = 0; k < n; k++) {
				dtheta1 = weight_theta[k];
				theta1 = theta[k];
				for (int l = 0; l < n; l++) {
					dtheta2 = weight_theta[l];
					theta2 = theta[l];
					for (int m = 0; m < n; m++) {
						dphi1 = weight_phi[m];
						phi1 = phi[m];
						for (int o = 0; o < n; o++) {
							dphi2 = weight_phi[o];
							phi2 = phi[o];
							integral += dr1*dr2*dtheta1*dtheta2*dphi1*dphi2 * func_polar(r1, theta1, phi1, r2, theta2, phi2);
						}
					}
				}
			}
		}
	}
	double exact = 5*pi*pi/(16.0*16.0);
    cout << "N          = " << n << endl;
    cout << "Numeric    = " << integral << endl;
    cout << "Exact      = " << exact << endl;
    cout << "Rel. error = " << (abs(exact-integral))/exact << endl;    
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    cout << "Time used for brute force polar: " << duration*1e-6 << " sec." << endl; 
}
void functions::monte_carlo_cartesian(int n, double a, double b){
    double const pi = 3.14159265359;
    double x1; double y1; double z1; double x2; double y2; double z2;
    double MCint, MCintsqr2, fx, Variance; 
    MCint = MCintsqr2=0.;
    double invers_period = 1./RAND_MAX;
    srand(time(NULL));  
    auto t1 = std::chrono::high_resolution_clock::now();
    for ( int i = 1;  i <= n; i++){
        x1 =  a + (b-a)*double(rand())*invers_period;
        x2 =  a + (b-a)*double(rand())*invers_period;
        y1 =  a + (b-a)*double(rand())*invers_period;
        y2 =  a + (b-a)*double(rand())*invers_period;
        z1 =  a + (b-a)*double(rand())*invers_period;
        z2 =  a + (b-a)*double(rand())*invers_period;
        fx = func_cartesian(x1, y1, z1, x2, y2, z2);
        MCint += fx;
        MCintsqr2 += fx*fx;
    }
    MCint = MCint/((double) n );
    MCintsqr2 = MCintsqr2/((double) n );
    double variance = MCintsqr2 - MCint*MCint;
    double final_MCint = MCint*pow((b-a),6);
    double exact = 5*pi*pi/(16.0*16.0);
    cout << "Variance   = " << variance << endl;
    cout << "Numeric    = " << final_MCint << endl;
    cout << "Exact      = " << exact << endl;
    cout << "Rel. error = " << (abs(exact-final_MCint))/exact << endl;
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    cout << "Time used for Monte Carlo Cartesian: " << duration*1e-6 << " sec." << endl;
}
void functions::monte_carlo_polar(int n){
    double const pi = 3.14159265359;
    double r1; double theta1; double phi1; double r2; double theta2; double phi2;
    double MCint, MCintsqr2, fx, Variance; 
    MCint = MCintsqr2 = 0.;
    double invers_period = 1./RAND_MAX;
    auto t1 = std::chrono::high_resolution_clock::now();
    srand(time(NULL));  
    for ( int i = 1;  i <= n; i++){
        r1 =  -log(1 - double(rand())*invers_period);
        r2 =  -log(1 - double(rand())*invers_period);
        theta1 =  pi*double(rand())*invers_period;
        theta2 =  pi*double(rand())*invers_period;
        phi1 =  2*pi*double(rand())*invers_period;
        phi2 =  2*pi*double(rand())*invers_period;
        fx = func_polar(r1, theta1, phi1, r2, theta2, phi2);
        MCint += fx;
        MCintsqr2 += fx*fx;
    }
    MCint = MCint/((double) n );
    MCintsqr2 = MCintsqr2/((double) n );
    double variance = MCintsqr2 - MCint*MCint;
    double powpow = 4*pow(pi, 4);
    double final_MCint = MCint*powpow;
    double exact = 5*pi*pi/(16.0*16.0);
    cout << "Variance   = " << variance << endl;
    cout << "Numeric    = " << final_MCint << endl;
    cout << "Exact      = " << exact << endl;
    cout << "Rel. error = " << (abs(exact-final_MCint))/exact << endl;
    auto t2 = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    cout << "Time used for Monte Carlo polar: " << duration*1e-6 << " sec." << endl;   

}
void functions::MPI_monte_carlo_cartesian(int n, double a, double b){
    int local_n, numprocs, my_rank; 
	double  h, local_a, local_b, total_sum;   
	double  time_start, time_end, total_time;
    int nargs = 1;
    char *args[] = {};

	MPI_Init (NULL,NULL);
	MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
	time_start = MPI_Wtime();
	total_sum = 0.0;
    
    double const pi = 3.14159265359;
    double x1; double y1; double z1; double x2; double y2; double z2;
    double MCint, fx; 
    MCint = 0.;
    double invers_period = 1./RAND_MAX;
    srand(time(NULL)*(my_rank+1));  
    for ( int i = 1;  i <= n; i++){
        x1 =  a + (b-a)*double(rand())*invers_period;
        x2 =  a + (b-a)*double(rand())*invers_period;
        y1 =  a + (b-a)*double(rand())*invers_period;
        y2 =  a + (b-a)*double(rand())*invers_period;
        z1 =  a + (b-a)*double(rand())*invers_period;
        z2 =  a + (b-a)*double(rand())*invers_period;
        fx = func_cartesian(x1, y1, z1, x2, y2, z2);
        MCint += fx;
    }
    double final_MCint = (MCint/((double) n ))*pow((b-a),6);
    double exact = 5*pi*pi/(16.0*16.0);

    MPI_Reduce(&final_MCint, &total_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	time_end = MPI_Wtime();
	total_time = time_end-time_start;
    //cout << "hei" << my_rank << endl;
	if ( my_rank == 0) {
        cout << "Monte-Carlo integration, cartesian integrand" << endl;
        cout << "Numeric    = " << total_sum/4 << endl;
        cout << "Exact      = " << exact << endl;   
        cout << "Rel. error = " << (abs(exact-final_MCint))/exact << endl;
		cout << "Time = " <<  total_time  << " on number of processors: "  << numprocs  << endl;
	}
	MPI_Finalize();  
}
void functions::MPI_monte_carlo_polar(int n){
    int local_n, numprocs, my_rank; 
	double  h, total_sum;   
	double  time_start, time_end, total_time;
    int nargs = 1;
    char *args[] = {};

	MPI_Init (NULL,NULL);
	MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
	time_start = MPI_Wtime();
	total_sum = 0.0;
    
    double const pi = 3.14159265359;
	double r1; double theta1; double phi1; double r2; double theta2; double phi2;
    double MCint, fx, Variance; 
    MCint = 0.;
    double invers_period = 1./RAND_MAX;
	srand(time(NULL)*(my_rank+1));  
    for ( int i = 1;  i <= n; i++){
        r1 =  -log(1 - double(rand())*invers_period);
        r2 =  -log(1 - double(rand())*invers_period);
        theta1 =  pi*double(rand())*invers_period;
        theta2 =  pi*double(rand())*invers_period;
        phi1 =  2*pi*double(rand())*invers_period;
        phi2 =  2*pi*double(rand())*invers_period;
        fx = func_polar(r1, theta1, phi1, r2, theta2, phi2);
        MCint += fx;
    }
    double powpow = 4*pow(pi, 4);
    double final_MCint = MCint*powpow;
    double exact = 5*pi*pi/(16.0*16.0);

    MPI_Reduce(&final_MCint, &total_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	time_end = MPI_Wtime();
	total_time = time_end-time_start;
    //cout << "hei" << my_rank << endl;
	if ( my_rank == 0) {
        cout << "Monte-Carlo integration, polar integrand" << endl;
        cout << "Numeric    = " << total_sum/4 << endl;
        cout << "Exact      = " << exact << endl;   
        cout << "Rel. error = " << (abs(exact-final_MCint))/exact << endl;
		cout << "Time = " <<  total_time  << " on number of processors: "  << numprocs  << endl;
	}
	MPI_Finalize();  

}
void functions::tester_bestfunc(){
    int n = 1e7;
    double a = -3;
    double b = 3;
    double const pi = 3.14159265359;
    double r1; double theta1; double phi1; double r2; double theta2; double phi2;
    double MCint, fx; 
    MCint = 0.;
    double invers_period = 1./RAND_MAX;
    srand(time(NULL));  
    for ( int i = 1;  i <= n; i++){
        r1 =  -log(1 - double(rand())*invers_period);
        r2 =  -log(1 - double(rand())*invers_period);
        theta1 =  pi*double(rand())*invers_period;
        theta2 =  pi*double(rand())*invers_period;
        phi1 =  2*pi*double(rand())*invers_period;
        phi2 =  2*pi*double(rand())*invers_period;
        fx = func_polar(r1, theta1, phi1, r2, theta2, phi2);
        MCint += fx;
    }
    MCint = MCint/((double) n );
    double powpow = 4*pow(pi, 4);
    double final_MCint = MCint*powpow;
    double exact = 5*pi*pi/(16.0*16.0);
    double tolerance = 1e-4;
    if (abs(exact - final_MCint)<tolerance){
        cout << "Test for checing the exact value of the integral is passed." << endl;
    }
    else{
        cout << "Something wierd happen, the calculated values and the exact values are not the same." << endl;
        cout << "Integral    = " << final_MCint << endl;
        cout << "Exact       = " << exact << endl;  
    }
}
void functions::tester_computed_exact_values(){
    double const pi = 3.14159265359;
    double exact = 5*pi*pi/(16.0*16.0);
    double value = 0.192765711;
    double tolerance = 1e-4;
    if (abs(exact - value)<tolerance){
        cout << "Test for checing the exact value of the integral is passed." << endl;
    }
    else{
        cout << "Something wierd happen, the calculated values and the exact values are not within the tolerance." << endl;
        cout << "Value    = " << value << endl;
        cout << "Exact    = " << exact << endl;
    }
}