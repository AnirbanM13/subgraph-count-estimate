# Load RSpectra package for fast partial eigen-decomposition
library(RSpectra)

# Compile and load the C++ Bethe Hessian implementation
Rcpp::sourceCpp("Desktop/r-file/bethe_hessian.cpp")

# Define Bethe Hessian spectral clustering function
bethe_hessian_clustering_adj <- function(A,
                                         K_opt = NULL,
                                         threads = 4,
                                         max_iter = 100) {
  
  N <- nrow(A)  # Number of nodes (assumes A is NxN adjacency matrix)
  
  # ---- degree & r ----
  
  deg <- rowSums(A)        # Compute node degrees
  r <- sqrt(mean(deg))     # Compute Bethe-Hessian regularization parameter:
  # r = sqrt(average degree)
  
  # ---- number of clusters ----
  
  # If K_opt not provided, default to 2 clusters
  K <- if (is.null(K_opt)) 2 else K_opt
  
  # ---- RSpectra operator (MUST accept 2 args) ----
  
  # Define matrix-vector multiplication operator
  # RSpectra does not need full matrix explicitly —
  # it only needs a function that computes H %*% x
  op <- function(x, args) {
    bethe_hessian_matvec(A, deg, r, x, threads)
    # Calls C++ function to compute:
    # H(r) x = ((r^2 - 1)I - rA + D) x
  }
  
  # ---- partial eigendecomposition ----
  
  # Compute K smallest algebraic eigenvalues/eigenvectors
  # of the Bethe Hessian matrix using RSpectra
  eig <- eigs_sym(
    A = op,        # Matrix operator
    k = K,         # Number of eigenvectors to compute
    n = N,         # Matrix dimension
    which = "SA"   # Smallest algebraic eigenvalues
  )
  
  X <- eig$vectors  # Extract eigenvectors (NxK matrix)
  
  # ---- k-means ----
  
  # Perform k-means clustering on rows of eigenvector matrix
  km <- kmeans(X, centers = K, iter.max = max_iter, nstart = 10)
  # centers = K clusters
  # iter.max = maximum iterations
  # nstart = run kmeans 10 times for stability
  
  km$cluster  # Return cluster labels
}