#include "operations.hpp"
#include "timer.hpp"
#include <omp.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <numeric>
#include <algorithm>
#include <cmath>

int main(){

int max_threads = omp_get_max_threads();

std::cout << "Maximum number of threads: " << max_threads << std::endl;
std::cout << "STRONG CONVERGENCE: " << std::endl;
std::cout << "--------------------" << std::endl;

// Strong scaling experiments: 
// Problem size constant, number of nodes/processors/cores increased

const long int n = 300000; // the size of the input array
double* x = new double[n]; // allocate the input array
double value = 1.0; // the value to fill the array with

// run the function with different numbers of threads
// strong covergence for init()
for (int num_threads = 1; num_threads <= max_threads; num_threads++)
{
// set the number of threads
omp_set_num_threads(num_threads);

std::cout << "Time taken by init() for problem size n = " << n << " and threads/cores p = " << num_threads  << std::endl;
Timer timer("init-strong");
init(n, x, value);

Timer::summarize(std::cout);
timer.~Timer();
}

double* y = new double[n]; // allocate the y array

for (int i=0; i<n; i++)
{
x[i] = double(i+1);
y[i] = 1.0/double(i+1);
}

// run the function with different numbers of threads
// strong convergence for dot()
for (int num_threads = 1; num_threads <= max_threads; num_threads++)
{
// set the number of threads
omp_set_num_threads(num_threads);

std::cout << "Time taken by dot() for problem size n = " << n << " and threads/cores p = " << num_threads  << std::endl;
Timer timer("dot-strong");
dot(n, x, y);
Timer::summarize(std::cout);
    
timer.~Timer();
}

double val1 = 1.0;
double val2 = 2.0;
init(n,x,val1);
init(n,y,val2);
double a = 3.0;
double b = 4.0;

// run the function with different numbers of threads
// strong convergence for axpby()
for (int num_threads = 1; num_threads <= max_threads; num_threads++)
{
// set the number of threads
omp_set_num_threads(num_threads);

std::cout << "Time taken by axpby() for problem size n = " << n << " and threads/cores p = " << num_threads  << std::endl;
Timer timer("axpby-strong");
axpby(n,a,x,b,y);
Timer::summarize(std::cout);
    
timer.~Timer();
}

// creation of appropriate 7-point stencil

int Nx = int(n/3);
int Ny = int(n/3);
int Nz = int(n/3);
const double dx = 1.0/(Nx-1);
const double dy = 1.0/(Ny-1);
const double dz = 1.0/(Nz-1);

stencil3d S;

S.nx=Nx; S.ny=Ny; S.nz=Nz;
S.value_c = 2.0/(dx*dx)+ 2.0/(dy*dy) + 2.0/(dz*dz);
S.value_n = -1.0/(dx*dx);
S.value_e = -1.0/(dx*dx);
S.value_s = -1.0/(dy*dy);
S.value_w = -1.0/(dy*dy);
S.value_b = -1.0/(dz*dz);
S.value_t = -1.0/(dz*dz);

// run the function with different numbers of threads
// strong convergence for apply_stencil3d()
for (int num_threads = 1; num_threads <= max_threads; num_threads++)
{
    
// Initialize u vector
double* u=new double[n];
for (int i=0; i<n; i++) u[i]=double(i);

// Create vector v where the solutions of Ae = v will be stored (for each iteration i)
double* v=new double[n*n];
//init(n,v,0.0);
omp_set_num_threads(num_threads); // set number of threads

std::cout << "Time taken by apply_stencil3d() for problem size n = " << n << " and threads/cores p = " << num_threads  << std::endl;
Timer timer("apply_stencil3d");
apply_stencil3d(&S, u, v); 
Timer::summarize(std::cout);

timer.~Timer();

delete[] u; // free memory
delete[] v; // free memory
}


// Weak scaling experiments:
// Problem size increased propertionally with the number of nodes/processors/cores

std::cout << "WEAK CONVERGENCE: " << std::endl;
std::cout << "--------------------" << std::endl;

const int min_n = 300000;    // min problem size
const int max_n = 3000000;   // max problem size
const int min_p = 1;        // min number of processors
const int max_p = max_threads;       // max number of processors

// std::ofstream out-init-weak("init-weak.txt");

// weak convergence for init()
for(int p = min_p; p <= max_p; p++) // loop over number of cores
    {
        for(int n = min_n; n <= max_n; n *= 10) // loop over problem size
        {
            double* x = new double[n]; // array to initialize
            double value = 1.0; // value to initialize array with

            omp_set_num_threads(p); // set number of threads
            
            std::cout << "Time taken by init() for problem size n = " << n << " and threads/cores p = " << p  << std::endl;
            Timer timer("init");
            init(n, x, value);
            // Timer::summarize(out-init-weak) // end timer
            Timer::summarize(std::cout);
            
            timer.~Timer();
            delete[] x; // free memory
        }
    }

// weak convergence for dot()    
for(int p = min_p; p <= max_p; p++) // loop over number of cores
    {
        for(int n = min_n; n <= max_n; n *= 10) // loop over problem size
        {
            double* x = new double[n]; // x array
            double* y = new double[n]; // y array

            omp_set_num_threads(p); // set number of threads
            
            std::cout << "Time taken by dot() for problem size n = " << n << " and threads/cores p = " << p  << std::endl;
            Timer timer("dot");
            dot(n, x, y);
            // Timer::summarize(out-init-weak) // end timer
            Timer::summarize(std::cout);
            
            timer.~Timer();
            delete[] x; // free memory
            delete[] y; // free memory
        }
    }

// weak convergence for axpby()    
for(int p = min_p; p <= max_p; p++) // loop over number of cores
    {
        for(int n = min_n; n <= max_n; n *= 10) // loop over problem size
        {
            double* x = new double[n]; // x array
            double* y = new double[n]; // y array
            double val1 = 1.0;
            double val2 = 2.0;
            init(n,x,val1);
            init(n,y,val2);
            double a = 3.0;
            double b = 4.0;

            omp_set_num_threads(p); // set number of threads
            
            std::cout << "Time taken by axpby() for problem size n = " << n << " and threads/cores p = " << p  << std::endl;
            Timer timer("axpby");
            axpby(n,a,x,b,y);
            // Timer::summarize(out-init-weak) // end timer
            Timer::summarize(std::cout);
            
            timer.~Timer();
            delete[] x; // free memory
            delete[] y; // free memory
        }
    }

// weak convergence for apply_stencil3d()    
for(int p = min_p; p <= max_p; p++) // loop over number of cores
    {
        for(int n = min_n; n <= max_n; n *= 10) // loop over problem size
        {
            
            int Nx = int(n/3);
            int Ny = int(n/3);
            int Nz = int(n/3);
            const double dx = 1.0/(Nx-1);
            const double dy = 1.0/(Ny-1);
            const double dz = 1.0/(Nz-1);
            
            stencil3d T;

            T.nx=Nx; T.ny=Ny; T.nz=Nz;
            T.value_c = 2.0/(dx*dx)+ 2.0/(dy*dy) + 2.0/(dz*dz);
            T.value_n = -1.0/(dx*dx);
            T.value_e = -1.0/(dx*dx);
            T.value_s = -1.0/(dy*dy);
            T.value_w = -1.0/(dy*dy);
            T.value_b = -1.0/(dz*dz);
            T.value_t = -1.0/(dz*dz);
  
            // Initialize u vector
            double* u=new double[n];
            for (int i=0; i<n; i++) u[i]=double(i);
  
            // Create vector v where the solutions of Ae = v will be stored (for each iteration i)
            double* v=new double[n*n];
            for (int i=0; i<n; i++) v[i]=0.0;

            omp_set_num_threads(p); // set number of threads
            
            std::cout << "Time taken by apply_stencil3d() for problem size n = " << n << " and threads/cores p = " << p  << std::endl;
            Timer timer("apply_stencil3d");
            apply_stencil3d(&T, u, v);
            Timer::summarize(std::cout);
            
            timer.~Timer();
            delete[] u; // free memory
            delete[] v; // free memory
        }
    }

return 0;
}