#include "gtest_mpi.hpp"

#include "cg_solver.hpp"

#include <iostream>

#include <vector>

#include <cmath>

#include <limits>

TEST(cg_solver, cg_solver){
    
  // TEST 1: Expected versus Numerical solution
  
  // Symmetric positive definite matrix, random numbers
  
  const int Nx=3, Ny=3, Nz=3;
  const int N=Nx*Ny*Nz;

  stencil3d S;

  S.nx=Nx; S.ny=Ny; S.nz=Nz;
  S.value_c = 1;
  S.value_n = 33;
  S.value_e = 5;
  S.value_s = 33;
  S.value_w = 5;
  S.value_b = 7;
  S.value_t = 7;
  
  // Set x_vector and initializing of b_vector
  
  double* x_vector = new double[N];
  for (int i=0; i<N; i++) x_vector[i]=i+1;
    
  double* b_vector = new double[N];
  init(N, b_vector, 0.0);
  
  apply_stencil3d(&S, x_vector, b_vector);
  
  // We have calculated Ax = b with our symmetric A and a vector x
  // We want to find the x back through applying the cg_solver() on A and b
  
  double* x_cg = new double[N];
  init(N, x_cg, 0.0);
  
  double tol = std::sqrt(std::numeric_limits<double>::epsilon());;
  int iter_max = 500;
  
  double res_norm;
  int iter_num;

  cg_solver(&S, N, x_cg, b_vector, 
  tol, iter_max, &res_norm, &iter_num);

  // the difference between the chosen x_vector and the obtained x_cg
  
  double err=0.0;
    for (int ii=0; ii<N; ii++) err = std::max(err, std::abs(x_vector[ii]-x_cg[ii]));
  
  EXPECT_NEAR(1.0+err, 1.0, 1e-10);
  
  // TEST 2: Correct number of iterations
  
  ASSERT_GT(N, iter_num); // the number of iterations should be equal or smaller than the matrix dimensions
  
  // 1D laplace problem
  
  /*
  const int nx=4, ny=0, nz=0;
  const int n=nx;

  stencil3d T;

  T.nx=nx; T.ny=ny; T.nz=nz;
  T.value_c = 2;
  T.value_n = 0;
  T.value_e = -1;
  T.value_s = 0;
  T.value_w = -1;
  T.value_b = 0;
  T.value_t = 0;
  
  double* b = new double[n];
  for (int j=0; j<n; j++) b_vector[j]=j+1;
  
  double* x = new double[n];
  init(n, x, 0.0);

  double tl = std::sqrt(std::numeric_limits<double>::epsilon());;
  int i_max = 20;
  
  double r_norm;
  int i_num;

  cg_solver(&T, n, x, b, 
  tl, i_max, &r_norm, &i_num);
  
*/
  
}