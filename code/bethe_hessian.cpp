// Enable OpenMP support for parallel execution
// [[Rcpp::plugins(openmp)]]

// Tell Rcpp to link against RcppEigen
// [[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>   // Rcpp + Eigen linear algebra integration
#include <omp.h>         // OpenMP header for multithreading
#include <cmath>         // Standard math functions (e.g., sqrt)

using namespace Rcpp;    // Avoid writing Rcpp::
using namespace Eigen;   // Avoid writing Eigen::

/*
 Matrix-vector product for Bethe Hessian:
 
 Hx = (r^2 - 1)x + D x - r A x
 
 where:
 A = adjacency matrix
 D = diagonal degree matrix
 r = regularization parameter
 x = input vector
 y = output vector
 */

// Export function to R
// [[Rcpp::export]]
Eigen::VectorXd bethe_hessian_matvec(const Eigen::MatrixXd& A,
                                     const Eigen::VectorXd& deg,
                                     double r,
                                     const Eigen::VectorXd& x,
                                     int threads = 4) {
  
  // Number of nodes (matrix dimension)
  int N = A.rows();
  
  // Output vector of size N
  Eigen::VectorXd y(N);
  
  // Set number of OpenMP threads
  omp_set_num_threads(threads);
  
  // Parallel loop over rows
#pragma omp parallel for
  for (int i = 0; i < N; i++) {
    
    // Compute (A x)_i = dot product of row i of A with vector x
    double Ax_i = A.row(i).dot(x);
    
    // Compute Bethe Hessian product:
    // y_i = (r^2 - 1) * x_i
    //       + deg_i * x_i
    //       - r * (A x)_i
    y[i] = (r * r - 1.0) * x[i]
    + deg[i] * x[i]
    - r * Ax_i;
  }
  
  // Return resulting vector
  return y;
}