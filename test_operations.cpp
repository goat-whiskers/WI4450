#include "gtest_mpi.hpp"

#include "operations.hpp"

#include <iostream>

#include <vector>

#include <cmath>

// note: you may add any number of tests to verify
// your code behaves correctly, but do not change
// the existing tests.

TEST(stencil, bounds_check)
{
  stencil3d S;
  S.nx=5;
  S.ny=3;
  S.nz=2;
  EXPECT_THROW(S.index_c(-1,0,0), std::runtime_error);
  EXPECT_THROW(S.index_c(S.nx,0,0), std::runtime_error);
  EXPECT_THROW(S.index_c(0,-1,0), std::runtime_error);
  EXPECT_THROW(S.index_c(0,S.ny,0), std::runtime_error);
  EXPECT_THROW(S.index_c(0,0,-1), std::runtime_error);
  EXPECT_THROW(S.index_c(0,0,S.nz), std::runtime_error);
}

TEST(stencil, index_order_kji)
{
  stencil3d S;
  S.nx=50;
  S.ny=33;
  S.nz=21;

  int i=10, j=15, k=9;

  EXPECT_EQ(S.index_c(i,j,k), S.index_c(i-1,j,k)+1);
  EXPECT_EQ(S.index_c(i,j,k), S.index_c(i,j-1,k)+S.nx);
  EXPECT_EQ(S.index_c(i,j,k), S.index_c(i,j,k-1)+S.nx*S.ny);
}

TEST(operations, init)
{
    
  // create a vector of size 15, with numbers from 1 to 15 
  const int n=15;
  double x[n];
  for (int i=0; i<n; i++) x[i]=double(i+1);
  
  // initialize this vector with constant values
  double val=42.0;
  init(n, x, val);

  // the difference between the set value and the values inside the vector
  double err=0.0;
    for (int i=0; i<n; i++) err = std::max(err, std::abs(x[i]-val));

  // note: EXPECT_NEAR uses a tolerance relative to the size of the target,
  // near 0 this is very small, so we use an absolute test instead by 
  // comparing to 1 instead of 0.
  
  EXPECT_NEAR(1.0+err, 1.0, std::numeric_limits<double>::epsilon());
  
  // MY OWN UNIT TEST
  
  // initialize vector with init()
  const int N = 3;
  double vals = 8.0;
  double* z = new double[N];
  init(N, z, vals);
  
  // vector_size*element_values = sum of elements
  double result = vals*N;
  
  double sum = 0.0;
  
  for (int j = 0; j<N; j++){
      sum += z[j];
  }
  
  EXPECT_DOUBLE_EQ(sum, result);
}


TEST(operations, dot) {
  const int n=150;
  double x[n], y[n];

  for (int i=0; i<n; i++)
  {
    x[i] = double(i+1);
    y[i] = 1.0/double(i+1);
  }

  double res = dot(n, x, y);
  
  EXPECT_NEAR(res, (double)n, n*std::numeric_limits<double>::epsilon());
  
  // MY OWN UNIT TEST
  
  // initialize two vectors
    
    const int N = 17;
    double u[N], v[N];
    
    for (int i=0; i<N; i++){
        
    u[i] = double(i+1);
    v[i] = double(i+2);
    
    }
    
  // check symmetry, dot product property
  
  double dot1 = dot(N, u, v);
  double dot2 = dot(N, v, u);
  
  double err = std::abs(dot1 - dot2);
  
  EXPECT_NEAR(1.0+err, 1.0, std::numeric_limits<double>::epsilon());
  
}

TEST(operations, axpby) {

// MY OWN UNIT TEST

const int n=5;

double a = 2;
double b = 3;

double x[n];
double y[n];
double res[n];

for (int i=0; i<n; i++)
{
x[i] = double(i+1);
y[i] = 1.0 + double(i+1);
res[i] = 1.0 + double(i+1);
}

// distributive property test
// dot(n, x, axpby(n,a,x,b,y)) = a*dot(x,x) + b*dot(x,y)

axpby(n,a,x,b,res);

double lhs = dot(n, x, res);
double rhs = a*dot(n,x,x) + b*dot(n,x,y);

EXPECT_DOUBLE_EQ(rhs,lhs); // checks if the two values are (close to) the same

}

TEST(operations,stencil3d_symmetric)
{
//  const int nx=3, ny=4, nz=5;
  const int nx=2, ny=2, nz=2;
  const int n=nx*ny*nz;
double* e=new double[n];
for (int i=0; i<n; i++) e[i]=0.0;
  double* A=new double[n*n];

stencil3d S;

S.nx=nx; S.ny=ny; S.nz=nz;
S.value_c = 8;
S.value_n = 2;
S.value_e = 4;
S.value_s = 2;
S.value_w = 4;
S.value_b = 1;
S.value_t = 1;

for (int i=0; i<n; i++)
  {
    e[i]=1.0;
    if (i>0) e[i-1]=0.0;
    apply_stencil3d(&S, e, A+i*n);
  }

  int wrong_entries=0;
  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++)
    {
      if (A[i*n+j]!=A[j*n+i]) wrong_entries++;
    }
  EXPECT_EQ(0, wrong_entries);

  if (wrong_entries)
  {
    std::cout << "Your matrix (computed on a 2x2x2 grid by apply_stencil(I)) is ..."<<std::endl;
    for (int j=0; j<n; j++)
    {
      for (int i=0; i<n; i++)
      {
        std::cout << A[i*n+j] << " ";
      }
      std::cout << std::endl;
    }
  }
  delete [] e;
  delete [] A;
  
  
  // MY OWN UNIT TEST
  
  // Diagonal Dominance, only valid after the symmetry test passed
  // To make this universal, we would need to transpose the matrix somehow
  
  const int Nx=2, Ny=2, Nz=2;
  const int N=Nx*Ny*Nz;
  const double dx = 1.0/(Nx-1);
  const double dy = 1.0/(Ny-1);
  const double dz = 1.0/(Nz-1);

  // Standard 3D Laplace discretization stencil
  // The discretization matrix is symmetric, 
  // so the rows and columns with the same indices contain the same elements in the same order

  stencil3d T;

  T.nx=nx; T.ny=ny; T.nz=nz;
  T.value_c = 2.0/(dx*dx)+ 2.0/(dy*dy) + 2.0/(dz*dz);
  T.value_n = -1.0/(dx*dx);
  T.value_e = -1.0/(dx*dx);
  T.value_s = -1.0/(dy*dy);
  T.value_w = -1.0/(dy*dy);
  T.value_b = -1.0/(dz*dz);
  T.value_t = -1.0/(dz*dz);
  
  // Initialize unit vector ee
  // e will have 1.0 at index i and 0.0 everywhere else (for each iteration i)
  
  double* ee=new double[n];
  for (int i=0; i<N; i++) ee[i]=0.0;
  
  // Create vector v where the solutions of Ae = v will be stored (for each iteration i)
  
  double* v=new double[N*N];
  for (int i=0; i<N; i++) v[i]=0.0;
  
  double diag_element = 0.0; // where the main diagonal element will be stored
  double row_sum = 0.0; // where the sum of all off-diagonal in absolute values will be stored

  // for diagonal dominance we need |a_ii| > sum{|a_ij|}


  for (int h=0; h<N; h++)       // the number of "rows" inside the matrix 
  {
    ee[h]=1.0;
    if (h>0) ee[h-1]=0.0;
    apply_stencil3d(&T, ee, v);
    
    
    for (int i = 0; i<N; i++){
        
        std::cout << v[i] << std::endl;
        
        if (i == h){            // h = column index, i = row index of main diagonal element
            diag_element += v[i];
            //std::cout << v[i] << ", main-diagonal" << std::endl;
        }
    
        else{
        row_sum += fabs(v[i]);
            //std::cout << v[i] << ", off-diagonal" << std::endl;
        }
        
    }
    
    std::cout << h << ": " << "diag_element = " << diag_element << std::endl;
    std::cout << h << ": " << "row_sum = " << row_sum << std::endl;
    ASSERT_GT(diag_element,row_sum); // checks that diag_element > row_sum for all "rows"
    
    diag_element = 0.0;
    row_sum = 0.0;
    init(N,v,0.0);
    
  }

  
  
}
