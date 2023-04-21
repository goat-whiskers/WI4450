#include "operations.hpp"
#include "timer.hpp"
#include <omp.h>
#include <iostream>
#include <numeric>
#include <algorithm>

void init(int n, double* x, double value)
{

#pragma omp parallel for
for(int i = 0; i<n; i++)
{
x[i] = value;
}
return;
}

double dot(int n, double const* x, double const* y)
{
double sum_i = 0;

#pragma omp parallel for reduction(+:sum_i)
for (int i=0; i<n; i++)
{
sum_i += x[i]*y[i];
}

return sum_i;
}

void axpby(int n, double a, double const* x, double b, double* y)
{

#pragma omp parallel for
for(int i=0; i<n; i++)
{
y[i] = a*x[i] + b*y[i];
}
return;
}

//! apply a 7-point stencil to a vector
void apply_stencil3d(stencil3d const* S,
    double const* u, double* v)
{

double val_c, val_n, val_e, val_s, val_w, val_b, val_t;
int nx, ny, nz, n;

val_c = S->value_c;
val_n = S->value_n;
val_e = S->value_e;
val_s = S->value_s;
val_w = S->value_w;
val_b = S->value_b;
val_t = S->value_t;

nx = S->nx;
ny = S->ny;
nz = S->nz;

n = nx*ny*nz;

double* A_vec = new double[n];  // declare and allocate memory for the matrix row
init(n, A_vec, 0.0);  // initialize array with all zeros

double dot_product = 0;

#pragma omp parallel 
#pragma omp for ordered // collapse(3) ordered

// change order of for-loop order

for(int k = 0; k<nz; k++)
{

for(int j = 0; j<ny; j++)
{

    for(int i = 0; i<nx; i++)
        {
/*
for(int i = 0; i<nx; i++)
{

for(int j = 0; j<ny; j++)
{

    for(int k = 0; k<nz; k++)
        {
*/
        
         // creation of a row of matrix A for grid point with coordinate (i,j,k)
         
         #pragma omp ordered
         A_vec[S->index_c(i,j,k)] = val_c;  // we have a center value for all the grid points

        // with boundary elimination
    
        if(i<(nx-1)){
             A_vec[S->index_e(i,j,k)] = val_e; // elimination of eastern boundary
            }

        if(i>0){
            A_vec[S->index_w(i,j,k)] = val_w; // western boundary
            }   

        if(j<(ny-1)){
            A_vec[S->index_n(i,j,k)] = val_n; // northern boundary    
            }

        if(j>0){
            A_vec[S->index_s(i,j,k)] = val_s; // southern boundary    
            }

        if(k<(nz-1)){
            A_vec[S->index_t(i,j,k)] = val_t; // top boundary
            }

        if(k>0){
            A_vec[S->index_b(i,j,k)] = val_b; // bottom boundary
            }

        dot_product = dot(n, A_vec, u); // dot product of A_vec with u
        
        v[S->index_c(i,j,k)] = dot_product; // update v atomically
        
        init(n, A_vec, 0.0); // empty A_vec for the next iteration
        
    }
}

}

delete [] A_vec;


return;
}

// apply a Jacobi preconditioner to a vector
void apply_jacobi(stencil3d const* S,
double const* u, double* v)
{

double val_c;
int nx, ny, nz, n;

val_c = S->value_c;

nx = S->nx;
ny = S->ny;
nz = S->nz;

n = nx*ny*nz;

double* A_diag = new double[n];  // declare and allocate memory for the matrix diagonal
init(n, A_diag, 1/val_c);        // initialize diagonal array

// #pragma omp parallel for
for (int i=0; i<n; i++)
{
v[i] = A_diag[i]*u[i];
}

delete [] A_diag;

return;
}
