#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <chrono>
#define NDEBUG

#include "functions.h"
#define ZERO 1E-10;

using namespace std;

int main(){
	functions funcs;
	int n = 1e5;
	double a, b;
	a = -3, b = 3;

	//funcs.loopyloop_cartesian(n, a, b);
	//funcs.loopyloop_polar(n);
	
	//funcs.monte_carlo_cartesian(n, a, b);
	//funcs.monte_carlo_polar(n);
	//funcs.MPI_monte_carlo_cartesian(n, a, b);
	//funcs.MPI_monte_carlo_polar(n); 
	//funcs.tester_bestfunc();
	funcs.tester_computed_exact_values();
}
