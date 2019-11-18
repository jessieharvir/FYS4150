/* 
   Program to solve the two-dimensional Ising model 
   with zero external field using MPI
   The coupling constant J = 1
   Boltzmann's constant = 1, temperature has thus dimension energy
   Metropolis sampling is used. Periodic boundary conditions.
   The code needs an output file on the command line and the variables mcs, nspins,
   initial temp, final temp and temp step.
*/
#include "mpi.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
using namespace  std;

// output file
ofstream ofile;

// inline function for periodic boundary conditions
inline int periodic(int i, int limit, int add) { 
  return (i+limit+add) % (limit);
}
// Function to initialise energy and magnetization
void initialize(int, int **, double&, double&);
// The metropolis algorithm 
void Metropolis(int, long&, int **, double&, double&, double *, int&);
// prints to file the results of the calculations  
void output(int, int, double, double *, int, int, double);
//  Matrix memory allocation
//  allocate space for a matrix
void  **matrix(int, int, int);
//  free space for  a matrix
void free_matrix(void **);
// ran2 for uniform deviates, initialize with negative seed.
double ran2(long *);
void initialize_random(int n_spins, int **spin_matrix, double& E, double& M);
void print_screen(int n_spins, int mcs, double temperature, double *average);
void ising_model(int argc, char* argv[]);
void ising_model_print(int argc, char* argv[], int n_spins, int mcs, double temp);
void ising_model_mcs(int argc, char* argv[], int n_spins, int mcs, double temp, string outfilename);
void ising_model_temp(int argc, char* argv[], int n_spins, int mcs, double initial_temp, double final_temp, double temp_step, string outfilename);
// Main program begins here

int main(int argc, char* argv[]){
    string outfilename = "/home/jessie/Dokumenter/skole/H19/FYS4150/Project4/output_T10.txt";
    int n_spins = 20, mcs = 100000;
    //double temp = 2.4;
    double temp = 1.0;
    //ising_model_mcs(argc, argv, n_spins, mcs, temp, outfilename);
    //string outfilename = "/home/jessie/Dokumenter/skole/H19/FYS4150/Project4/output_L100_10e5_random.txt";
    //double initial_temp = 2.0, final_temp = 2.4, temp_step = 0.005;
    //ising_model_temp(argc, argv, n_spins, mcs, initial_temp, final_temp, temp_step, outfilename);

    //ising_model_print(argc, argv, n_spins, 1e2, temp);
    ising_model_print(argc, argv, n_spins, 1e3, temp);
    //ising_model_print(argc, argv, n_spins, 1e4, temp);
    //ising_model_print(argc, argv, n_spins, 1e5, temp);
    //ising_model_print(argc, argv, n_spins, 1e6, temp);


        
}
void ising_model_mcs(int argc, char* argv[], int n_spins, int mcs, double temp, string outfilename){
    long idum;
    int **spin_matrix, my_rank, numprocs;
    double w[17], average[5], total_average[5], E, M;
    ofile.open(outfilename);
    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
     
    int no_intervalls = mcs/numprocs;
    int myloop_begin = my_rank*no_intervalls + 1;
    int myloop_end = (my_rank+1)*no_intervalls;

    spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int));
    idum = -1-my_rank;
    double  TimeStart, TimeEnd, TotalTime;
    TimeStart = MPI_Wtime();
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(15) << setprecision(8) << "cycles";                           // monte carlo cycles
    ofile << setw(15) << setprecision(8) << "counter";                          // number of accepted configs 
    ofile << setw(15) << setprecision(8) << "temperature";                      // temperature
    ofile << setw(15) << setprecision(8) << "E avg";                            // avg energy
    ofile << setw(15) << setprecision(8) << "Heat cap";                         // heat capasity
    ofile << setw(15) << setprecision(8) << "M avg";                            // avg magnetisation
    ofile << setw(15) << setprecision(8) << "Magn susp";                        // magnetic suspect
    ofile << setw(15) << setprecision(8) << "M avg abs";                        // magnetic moment abs
    ofile << setw(15) << setprecision(8) << "E variance";                       // variance energy
    ofile << setw(15) << setprecision(8) << "M variance";                       // variance magnetic
    ofile << setw(15) << setprecision(8) << "E" << endl;                        // energy
    if ( (my_rank == numprocs-1) &&( myloop_end < mcs) ) myloop_end = mcs;
    E = M = 0.;
    initialize(n_spins, spin_matrix, E, M);
    for( int de =-8; de <= 8; de++) w[de+8] = 0;
    for( int de =-8; de <= 8; de+=4) w[de+8] = exp(-de/temp);
    for( int i = 0; i < 5; i++) average[i] = 0.;
    for( int i = 0; i < 5; i++) total_average[i] = 0.;
    int counter = 0;
    for (int cycles = myloop_begin; cycles <= myloop_end; cycles++){
        Metropolis(n_spins, idum, spin_matrix, E, M, w, counter);
        average[0] += E;    average[1] += E*E;
        average[2] += M;    average[3] += M*M; average[4] += fabs(M);
        output(n_spins, mcs, temp, average, counter, cycles, E);
        /*
        double Evariance = (average[1]- E*E)/n_spins/n_spins;
        ofile << setw(15) << setprecision(8) << cycles;     // monte carlo cycles
        ofile << setw(15) << setprecision(8) << E/n_spins/n_spins;     // avg energy
        ofile << setw(15) << setprecision(8) << M/n_spins/n_spins;     // magnetic moment
        ofile << setw(15) << setprecision(8) << fabs(M)/n_spins/n_spins; // magnetic moment abs
        ofile << setw(15) << setprecision(8) << counter; // magnetic moment abs   
        ofile << setw(15) << setprecision(8) << Evariance << endl; // magnetic moment abs   
        */

    }
    for( int i =0; i < 5; i++){
    MPI_Reduce(&average[i], &total_average[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    free_matrix((void **) spin_matrix); // free memory
    ofile.close();  // close output file
    TimeEnd = MPI_Wtime();
    TotalTime = TimeEnd-TimeStart;
    if ( my_rank == 0) {
        cout << "Time = " <<  TotalTime  << " on number of processors: "  << numprocs  << endl;
    }
    // End MPI
    MPI_Finalize (); 
}

void ising_model_print(int argc, char* argv[], int n_spins, int mcs, double temp){
    long idum;
    ising_model_print(argc, argv, n_spins, 1e2, temp);
    int **spin_matrix, my_rank, numprocs;
    ising_model_print(argc, argv, n_spins, 1e2, temp);
    double w[17], average[5], total_average[5], E, M;
    MPI_Init (&argc, &argv);
    ising_model_print(argc, argv, n_spins, 1e2, temp);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    ising_model_print(argc, argv, n_spins, 1e2, temp);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
     
    ising_model_print(argc, argv, n_spins, 1e2, temp);
    //  Allocate memory for spin matrix
    ising_model_print(argc, argv, n_spins, 1e2, temp);
    spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int));
    ising_model_print(argc, argv, n_spins, 1e2, temp);
    // every node has its own seed for the random numbers, this is important else
    // if one starts with the same seed, one ends with the same random numbers
    idum = -1-my_rank;  // random starting point
    // Start Monte Carlo sampling by looping over mcs first
    double  TimeStart, TimeEnd, TotalTime;
    TimeStart = MPI_Wtime();
        //    initialise energy and magnetization 
        E = M = 0.;
        // initialise array for expectation values
        initialize(n_spins, spin_matrix, E, M);
        // setup array for possible energy changes
        for( int de =-8; de <= 8; de++) w[de+8] = 0;
        for( int de =-8; de <= 8; de+=4) w[de+8] = exp(-de/temp);
        for( int i = 0; i < 5; i++) average[i] = 0.;
        for( int i = 0; i < 5; i++) total_average[i] = 0.;
        // start Monte Carlo computation
        int counter = 0;
        for (int cycles = mcs; cycles <= mcs; cycles++){
            Metropolis(n_spins, idum, spin_matrix, E, M, w, counter);
            // update expectation values  for local node
            average[0] += E;    average[1] += E*E;
            average[2] += M;    average[3] += M*M; average[4] += fabs(M);
        }
        // Find total average
        for( int i =0; i < 5; i++){
        MPI_Reduce(&average[i], &total_average[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }
        // print results
        if ( my_rank == 0) {
            print_screen(n_spins, mcs, temp, average);
        }
    free_matrix((void **) spin_matrix); // free memory
    TimeEnd = MPI_Wtime();
    TotalTime = TimeEnd-TimeStart;
    if ( my_rank == 0) {
        cout << "Time = " <<  TotalTime  << " on number of processors: "  << numprocs  << endl;
    }
    // End MPI
    MPI_Finalize (); 
}
/*
void ising_model_print(int argc, char* argv[], int n_spins, int mcs_start, int mcs_end, double temp, double mcs_step){
    long idum;
    int **spin_matrix, my_rank, numprocs;
    double w[17], average[5], total_average[5], E, M;
    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
     
    int no_intervalls = mcs_start/numprocs;
    int myloop_begin = my_rank*no_intervalls + 1;
    int myloop_end = (my_rank+1)*no_intervalls;
    // MPI initializations
    MPI_Bcast (&n_spins, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&mcs_start, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&mcs_end, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&mcs_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //  Allocate memory for spin matrix
    spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int));
    // every node has its own seed for the random numbers, this is important else
    // if one starts with the same seed, one ends with the same random numbers
    idum = -1-my_rank;  // random starting point
    // Start Monte Carlo sampling by looping over mcs first
    double  TimeStart, TimeEnd, TotalTime;
    TimeStart = MPI_Wtime();
    for (int mcs = mcs_start; mcs <= mcs_end;  mcs+=mcs_step){
        if ( (my_rank == numprocs-1) &&( myloop_end < mcs_end) ) myloop_end = mcs;
        //    initialise energy and magnetization 
        E = M = 0.;
        // initialise array for expectation values
        initialize(n_spins, spin_matrix, E, M);
        // setup array for possible energy changes
        for( int de =-8; de <= 8; de++) w[de+8] = 0;
        for( int de =-8; de <= 8; de+=4) w[de+8] = exp(-de/temp);
        for( int i = 0; i < 5; i++) average[i] = 0.;
        for( int i = 0; i < 5; i++) total_average[i] = 0.;
        // start Monte Carlo computation
        int counter = 0;
        for (int cycles = myloop_begin; cycles <= myloop_end; cycles++){
            Metropolis(n_spins, idum, spin_matrix, E, M, w, counter);
            // update expectation values  for local node
            average[0] += E;    average[1] += E*E;
            average[2] += M;    average[3] += M*M; average[4] += fabs(M);
        }
        // Find total average
        for( int i =0; i < 5; i++){
        MPI_Reduce(&average[i], &total_average[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }
        // print results
        if ( my_rank == 0) {
            print_screen(n_spins, mcs, temp, average);
        }
    }
    free_matrix((void **) spin_matrix); // free memory
    TimeEnd = MPI_Wtime();
    TotalTime = TimeEnd-TimeStart;
    if ( my_rank == 0) {
        cout << "Time = " <<  TotalTime  << " on number of processors: "  << numprocs  << endl;
    }
    // End MPI
    MPI_Finalize (); 
}*/


void ising_model_temp(int argc, char* argv[], int n_spins, int mcs, double initial_temp, double final_temp, double temp_step, string outfilename){
    long idum;
    int **spin_matrix, my_rank, numprocs;
    double w[17], average[5], total_average[5], E, M;
    ofile.open(outfilename);
    // MPI initializations
    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    int no_intervalls = mcs/numprocs;
    int myloop_begin = my_rank*no_intervalls + 1;
    int myloop_end = (my_rank+1)*no_intervalls;
    if ( (my_rank == numprocs-1) &&( myloop_end < mcs) ) myloop_end = mcs;

    // broadcast to all nodes common variables
    MPI_Bcast (&n_spins, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&initial_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&final_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&temp_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //  Allocate memory for spin matrix
    spin_matrix = (int**) matrix(n_spins, n_spins, sizeof(int));
    // every node has its own seed for the random numbers, this is important else
    // if one starts with the same seed, one ends with the same random numbers
    idum = -1-my_rank;  // random starting point
    // Start Monte Carlo sampling by looping over T first

    double  TimeStart, TimeEnd, TotalTime;
    TimeStart = MPI_Wtime();
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(15) << setprecision(8) << "cycles";                           // monte carlo cycles
    ofile << setw(15) << setprecision(8) << "counter";                          // number of accepted configs 
    ofile << setw(15) << setprecision(8) << "temperature";                      // temperature
    ofile << setw(15) << setprecision(8) << "E avg";                            // avg energy
    ofile << setw(15) << setprecision(8) << "Heat cap";                         // heat capasity
    ofile << setw(15) << setprecision(8) << "M avg";                            // avg magnetisation
    ofile << setw(15) << setprecision(8) << "Magn susp";                        // magnetic suspect
    ofile << setw(15) << setprecision(8) << "M avg abs";                        // magnetic moment abs
    ofile << setw(15) << setprecision(8) << "E variance";                       // variance energy
    ofile << setw(15) << setprecision(8) << "M variance";                       // variance magnetic
    ofile << setw(15) << setprecision(8) << "E" << endl;                        // energy
    for ( double temp = initial_temp; temp <= final_temp; temp+=temp_step){
        //    initialise energy and magnetization 
        E = M = 0.;
        // initialise array for expectation values
        initialize_random(n_spins, spin_matrix, E, M);
        // setup array for possible energy changes
        for( int de =-8; de <= 8; de++) w[de+8] = 0;
        for( int de =-8; de <= 8; de+=4) w[de+8] = exp(-de/temp);
        for( int i = 0; i < 5; i++) average[i] = 0.;
        for( int i = 0; i < 5; i++) total_average[i] = 0.;
        // start Monte Carlo computation
        int counter = 0;
        for (int cycles = myloop_begin; cycles <= myloop_end; cycles++){
            Metropolis(n_spins, idum, spin_matrix, E, M, w, counter);
            // update expectation values  for local node
            average[0] += E;    average[1] += E*E;
            average[2] += M;    average[3] += M*M; average[4] += fabs(M);
        }
        // Find total average
        for( int i =0; i < 5; i++){
            MPI_Reduce(&average[i], &total_average[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        }
        // print results
        if ( my_rank == 0) {
            int cycles = 1;
            output(n_spins, mcs, temp, total_average, counter, cycles, E);
        }
    }
    free_matrix((void **) spin_matrix); // free memory
    ofile.close();  // close output file
    TimeEnd = MPI_Wtime();
    TotalTime = TimeEnd-TimeStart;
    if ( my_rank == 0) {
        cout << "Time = " <<  TotalTime  << " on number of processors: "  << numprocs  << endl;
    }

    // End MPI
    MPI_Finalize (); 
}

// function to initialise energy, spin matrix and magnetization
void initialize(int n_spins, int **spin_matrix, double& E, double& M){
    // setup spin matrix and intial magnetization
    for(int y =0; y < n_spins; y++) {
        for (int x= 0; x < n_spins; x++){
        spin_matrix[y][x] = 1; // spin orientation for the ground state
        M +=  (double) spin_matrix[y][x];
        }
    }
    // setup initial energy
    for(int y =0; y < n_spins; y++) {
        for (int x= 0; x < n_spins; x++){
        E -=  (double) spin_matrix[y][x]*
        (spin_matrix[periodic(y,n_spins,-1)][x] +
        spin_matrix[y][periodic(x,n_spins,-1)]);
        }
    }
}// end function initialise

void initialize_random(int n_spins, int **spin_matrix, double& E, double& M){
    // setup spin matrix and intial magnetization
    for(int y =0; y < n_spins; y++) {
        for (int x= 0; x < n_spins; x++){
            int random = rand() %10;
            if (random < 5){
                spin_matrix[y][x] = -1; // spin orientation for the ground state
            }
            else{
                spin_matrix[y][x] = 1; // spin orientation for the ground state
            }
            M +=  (double) spin_matrix[y][x];
        }
    }
    // setup initial energy
    for(int y =0; y < n_spins; y++) {
        for (int x= 0; x < n_spins; x++){
            E -=  (double) spin_matrix[y][x]*
            (spin_matrix[periodic(y,n_spins,-1)][x] +
            spin_matrix[y][periodic(x,n_spins,-1)]);
        }
    }
}// end function initialise random

void Metropolis(int n_spins, long& idum, int **spin_matrix, double& E, double&M, double *w, int &counter){
    // loop over all spins
    for(int y =0; y < n_spins; y++) {
        for (int x= 0; x < n_spins; x++){
        int ix = (int) (ran2(&idum)*(double)n_spins);
        int iy = (int) (ran2(&idum)*(double)n_spins);
        int deltaE =  2*spin_matrix[iy][ix]*
        (spin_matrix[iy][periodic(ix,n_spins,-1)]+
        spin_matrix[periodic(iy,n_spins,-1)][ix] +
        spin_matrix[iy][periodic(ix,n_spins,1)] +
        spin_matrix[periodic(iy,n_spins,1)][ix]); 
        if ( ran2(&idum) <= w[deltaE+8] ) {
            spin_matrix[iy][ix] *= -1;  // flip one spin and accept new spin config
            M += (double) 2*spin_matrix[iy][ix];
            E += (double) deltaE;
            counter += 1;
            }
        }
    }
} // end of Metropolis sampling over spins


void output(int n_spins, int mcs, double temp, double *average, int counter, int cycles, double E){
    double norm = 1./((double) (mcs));  // divided by total number of cycles   
    double Etotal_average = average[0]*norm;
    double E2total_average = average[1]*norm;
    double Mtotal_average = average[2]*norm;
    double M2total_average = average[3]*norm;
    double Mabstotal_average = average[4]*norm;
    // all expectation values are per spin, divide by 1/n_spins/n_spins
    double Evariance = (E2total_average- Etotal_average*Etotal_average)/n_spins/n_spins;
    double Mvariance = (M2total_average - Mtotal_average*Mtotal_average)/n_spins/n_spins;
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << setw(15) << setprecision(8) << cycles;                         // monte carlo cycles
    ofile << setw(15) << setprecision(8) << counter;                        // number of accepted configs 
    ofile << setw(15) << setprecision(8) << temp;                           // temperature
    ofile << setw(15) << setprecision(8) << Etotal_average/n_spins/n_spins; // avg energy
    ofile << setw(15) << setprecision(8) << Evariance/temp/temp;            // heat capasity
    ofile << setw(15) << setprecision(8) << Mtotal_average/n_spins/n_spins; // magnetic moment
    ofile << setw(15) << setprecision(8) << Mvariance/temp;                 // magnetic suspect
    ofile << setw(15) << setprecision(8) << Mabstotal_average/n_spins/n_spins; // magnetic moment abs
    ofile << setw(15) << setprecision(8) << Evariance;                      // variance energy
    ofile << setw(15) << setprecision(8) << Mvariance;                      // variance magnetic
    ofile << setw(15) << setprecision(8) << E << endl;                      // energy
} // end output function

void print_screen(int n_spins, int mcs, double temperature, double *average){
    double norm = 1/((double) (mcs)); // divided by total number of cycles
    double Eaverage = average[0]*norm;
    double E2average = average[1]*norm;
    double Maverage = average[2]*norm;
    double M2average = average[3]*norm;
    double Mabsaverage = average[4]*norm;
    // all expectation values are per spin, divide by 1/n_spins/n_spins
    double Evariance = (E2average- Eaverage*Eaverage)/n_spins/n_spins;
    double Mvariance = (M2average - Maverage*Maverage)/n_spins/n_spins;
    double M2variance = (M2average - Mabsaverage*Mabsaverage)/n_spins/n_spins;
    double heat_capacaty = Evariance/(temperature*temperature);
    double magn_suscept = Mvariance/temperature;
    cout << "Monte Carlo cycles: " << mcs << endl;
    cout << "Average energy: " << Eaverage/n_spins/n_spins << endl;
    //cout << "Energy variance: " << Evariance/temperature/temperature << endl;
    //cout << "Magnetic variance: "<< M2variance/temperature << endl;
    cout << "Average magnetic moment: " << Mabsaverage/n_spins/n_spins << endl;
    cout << "Heat capacity: " << heat_capacaty << endl;
    cout << "Magnetic susceptibility: " << magn_suscept << endl;
}