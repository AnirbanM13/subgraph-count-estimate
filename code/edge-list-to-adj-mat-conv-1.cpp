#include <Rcpp.h>              // Rcpp interface to connect C++ with R
#include <unordered_map>       // Hash map container for fast key → value lookup
using namespace Rcpp;          // Avoid writing Rcpp:: repeatedly

// Export function so it can be called directly from R
// [[Rcpp::export]]
IntegerMatrix edge_list_to_adj_unique(DataFrame edges, bool directed = false) {
  
  // Check that the input DataFrame has exactly 2 columns (from, to)
  if (edges.size() != 2) 
    stop("DataFrame must have exactly 2 columns (from, to)");
  
  // Extract first column as source nodes
  IntegerVector from = edges[0];
  
  // Extract second column as destination nodes
  IntegerVector to   = edges[1];
  
  // Number of edges
  int M = from.size();
  
  // ---------------------------------------------------------
  // 1. Map unique vertex labels → contiguous indices (0..N-1)
  // ---------------------------------------------------------
  
  // Create hash map: original vertex label → new index
  std::unordered_map<int,int> remap;
  
  // Reserve memory for efficiency (at most 2*M unique endpoints)
  remap.reserve(2 * M);
  
  int idx = 0;  // Counter for assigning new indices
  
  // Loop over all edges
  for (int k = 0; k < M; ++k) {
    
    // If "from" vertex not yet seen, assign new index
    if (remap.find(from[k]) == remap.end())
      remap[from[k]] = idx++;
    
    // If "to" vertex not yet seen, assign new index
    if (remap.find(to[k]) == remap.end())
      remap[to[k]] = idx++;
  }
  
  // Total number of unique vertices
  int N = idx;
  
  // ---------------------------------------------------------
  // 2. Build adjacency matrix (N x N)
  // ---------------------------------------------------------
  
  // Create N × N integer matrix initialized to 0
  IntegerMatrix adj(N, N);
  
  // Loop through each edge
  for (int k = 0; k < M; ++k) {
    
    // Get remapped index of source vertex
    int u = remap[from[k]];
    
    // Get remapped index of destination vertex
    int v = remap[to[k]];
    
    // Skip self-loops (u == v)
    if (u == v) continue;
    
    // Set edge u → v
    adj(u, v) = 1;
    
    // If graph is undirected, also set v → u
    if (!directed)
      adj(v, u) = 1;
  }
  
  // Return adjacency matrix to R
  return adj;
}