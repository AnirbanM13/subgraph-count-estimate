// Enable OpenMP support for parallel computation
// [[Rcpp::plugins(openmp)]]

#include <Rcpp.h>   // Rcpp interface to connect C++ with R
#include <omp.h>    // OpenMP library for multithreading

using namespace Rcpp;  // Avoid writing Rcpp:: repeatedly

// Export this function to R
// [[Rcpp::export]]
NumericMatrix PI_hat_cpp(IntegerVector a, NumericMatrix x, int threads = 1) {
  
  int N = x.nrow();  // Number of nodes (rows of adjacency matrix)
  
  // Ensure adjacency matrix is square
  if (x.ncol() != N) stop("Adjacency matrix must be square.");
  
  // Ensure cluster vector has correct length
  if (a.size() != N) stop("Length of a must equal dimension of x.");
  
  // Set number of OpenMP threads
  omp_set_num_threads(threads);
  
  // Determine number of clusters K = maximum label in a
  int K = *std::max_element(a.begin(), a.end());
  
  // Ensure cluster labels are positive
  if (K <= 0) stop("Cluster assignments must be positive integers.");
  
  // Convert cluster labels from R indexing (1..K) to C++ indexing (0..K-1)
  std::vector<int> cluster(N);  // Store zero-based cluster labels
  
  for (int i = 0; i < N; i++) {
    if (a[i] < 1) stop("Cluster labels must start from 1.");
    cluster[i] = a[i] - 1;  // Convert to zero-based indexing
  }
  
  // Count number of nodes in each cluster/block
  std::vector<long long> N_vec(K, 0);  // N_vec[r] = number of nodes in cluster r
  
  for (int i = 0; i < N; i++)
    N_vec[cluster[i]]++;  // Increment cluster count
  
  // Temporary storage for counting edges between blocks
  // Stored as flattened K x K matrix
  std::vector<long long> temp(K * K, 0);
  
  // Parallel loop over nodes
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < N; i++) {
    
    int ci = cluster[i];  // Cluster of node i
    
    // Only loop over upper triangle (j > i) to avoid double counting
    for (int j = i + 1; j < N; j++) {
      
      // Skip if no edge between i and j
      if (x(i, j) == 0) continue;
      
      int cj = cluster[j];  // Cluster of node j
      
      // Atomically increment block-to-block edge count
#pragma omp atomic
      temp[ci * K + cj] += 1;
      
      // Also increment symmetric entry (since graph is undirected)
#pragma omp atomic
      temp[cj * K + ci] += 1;
    }
  }
  
  // Create resulting K x K probability matrix
  NumericMatrix PI(K, K);
  
  // Compute block probabilities
  for (int r = 0; r < K; r++) {
    for (int s = r; s < K; s++) {
      
      // Off-diagonal blocks (r ≠ s)
      if (r != s) {
        
        // Ensure both blocks contain nodes
        if (N_vec[r] > 0 && N_vec[s] > 0) {
          
          // Estimate probability as:
          // (# edges between r and s) / (possible pairs)
          PI(r, s) = (double)temp[r * K + s] /
            ((double)N_vec[r] * (double)N_vec[s]);
          
          // Symmetric entry
          PI(s, r) = PI(r, s);
          
        } else {
          // If one block empty, probability = 0
          PI(r, s) = 0.0;
          PI(s, r) = 0.0;
        }
        
      } else {  // Diagonal block (within same cluster)
        
        // Need at least 2 nodes to form an edge
        if (N_vec[r] > 1) {
          
          // Within-block probability:
          // (# observed edges) / (N_r * (N_r - 1))
          // Note: since edges counted twice (symmetric),
          // denominator matches doubled counting.
          PI(r, r) = (double)temp[r * K + r] /
            ((double)N_vec[r] * (N_vec[r] - 1));
          
        } else {
          // If only one node, no possible edges
          PI(r, r) = 0.0;
        }
      }
    }
  }
  
  // Return estimated SBM probability matrix
  return PI;
}