// Enable OpenMP support
// [[Rcpp::plugins(openmp)]]

#include <Rcpp.h>          // Rcpp interface for R <-> C++
#include <random>          // Random number generators
#include <vector>          // std::vector container
#include <unordered_set>   // Hash set (used for unique edge sampling)
#include <cinttypes>       // Fixed-width integer types (uint64_t, etc.)
#include <omp.h>           // OpenMP parallelization
#include <algorithm>       // sort, min, max

using namespace Rcpp;      // Avoid writing Rcpp:: everywhere


/* =========================
 Utilities
 ========================= */

// Pack two 32-bit integers (u,v) into one 64-bit integer
// Used for fast uniqueness checking in hash set
static inline uint64_t pack_uv(uint32_t u, uint32_t v) {
  // Shift u into high 32 bits and OR with v
  return ((uint64_t)u << 32) | (uint64_t)v;
}


/* =========================
 Sort adjacency ONCE
 ========================= */

// Sort each adjacency list in parallel
void sort_adj(std::vector<std::vector<int>> &adj, int threads) {
  
  // Set number of OpenMP threads
  omp_set_num_threads(threads);
  
  // Parallel loop over all nodes
#pragma omp parallel for schedule(dynamic,256)
  for (int i = 0; i < (int)adj.size(); ++i)
    
    // Sort adjacency list of node i
    std::sort(adj[i].begin(), adj[i].end());
}


/* =========================
 Fast wedge + triangle count
 (assumes sorted adjacency)
 ========================= */

// Count wedges and triangles assuming adjacency lists are sorted
std::pair<long long,long long>
count_wedges_triangles_sorted(
  const std::vector<std::vector<int>> &adj,
  int threads) {
  
  int N = adj.size();   // Number of nodes
  
  long long wedges = 0, triangles = 0;  // Counters
  
  // ---- Count wedges ----
#pragma omp parallel for reduction(+:wedges)
  for (int i = 0; i < N; ++i) {
    
    // Degree of node i
    long long d = adj[i].size();
    
    // Add choose(d,2) if degree >= 2
    if (d >= 2) wedges += d * (d - 1) / 2;
  }
  
  // ---- Count triangles using sorted intersection ----
#pragma omp parallel for reduction(+:triangles) schedule(dynamic,50)
  for (int i = 0; i < N; ++i) {
    
    const auto &Ni = adj[i];  // Neighbors of i
    size_t di = Ni.size();    // Degree of i
    
    // Iterate neighbors j
    for (int j : Ni) {
      
      // Enforce ordering i < j
      if (j <= i) continue;
      
      const auto &Nj = adj[j];  // Neighbors of j
      
      // Degree pruning optimization:
      // If deg(i) > deg(j), skip to reduce comparisons
      if (di > Nj.size()) continue;
      
      // Two-pointer intersection
      size_t a = 0, b = 0;
      
      while (a < Ni.size() && b < Nj.size()) {
        
        if (Ni[a] == Nj[b]) {
          
          // Enforce ordering i < j < k
          if (Ni[a] > j) triangles++;
          
          a++; 
          b++;
          
        } else if (Ni[a] < Nj[b]) {
          a++;
        } else {
          b++;
        }
      }
    }
  }
  
  // Return (wedges, triangles)
  return {wedges, triangles};
}


/* =========================
 MAIN FUNCTION
 ========================= */

// Export to R
// [[Rcpp::export]]
List generate_sbm_edges_with_multiple_W_fast(
    IntegerVector alpha,      // Block assignments
    NumericMatrix Pi,         // Block probability matrix
    NumericVector probs,      // W probabilities
    int threads = 4,          // Number of threads
    int seed = -1) {          // RNG seed
  
  // Set OpenMP thread count
  omp_set_num_threads(threads);
  
  int N = alpha.size();   // Number of nodes
  int K = Pi.nrow();      // Number of blocks
  
  /* -------- group nodes -------- */
  
  // Create vector of node lists per block
  std::vector<std::vector<int>> groups(K);
  
  // Assign nodes to their block
  for (int i = 0; i < N; ++i)
    groups[alpha[i]-1].push_back(i);
  
  // Base random seed
  unsigned int base_seed =
    (seed < 0) ? std::random_device{}() : (unsigned int)seed;
  
  /* -------- generate SBM edges -------- */
  
  // Store generated edges
  std::vector<std::pair<int,int>> edges;
  edges.reserve(1000000);  // Pre-allocate memory
  
  // Full adjacency list
  std::vector<std::vector<int>> adj_all(N);
  
#pragma omp parallel
{
  // Thread-local RNG
  std::mt19937 gen(base_seed + 7919 * omp_get_thread_num());
  
  // Thread-local edge buffer
  std::vector<std::pair<int,int>> local_edges;
  local_edges.reserve(4096);
  
#pragma omp for schedule(dynamic)
  for (int a = 0; a < K; ++a) {
    for (int b = a; b < K; ++b) {
      
      auto &A = groups[a];
      auto &B = groups[b];
      
      // Skip empty blocks
      if (A.empty() || B.empty()) continue;
      
      double p = Pi(a,b);  // Edge probability
      if (p <= 0) continue;
      
      uint64_t na = A.size();
      uint64_t nb = B.size();
      
      // Total possible pairs
      uint64_t total =
        (a==b) ? na*(na-1)/2 : na*nb;
      
      // Sample number of edges M ~ Binomial(total, p)
      std::binomial_distribution<uint64_t> binom(total, p);
      uint64_t M = binom(gen);
      
      if (M == 0) continue;
      
      // Hash set to ensure uniqueness
      std::unordered_set<uint64_t> seen;
      seen.reserve(M*2 + 8);
      seen.max_load_factor(0.7f);
      
      std::uniform_int_distribution<uint32_t> da(0, na-1);
      std::uniform_int_distribution<uint32_t> db(0, nb-1);
      
      // Sample M unique edges
      while (seen.size() < M) {
        
        uint32_t u, v;
        
        if (a == b) {
          
          uint32_t i = da(gen);
          uint32_t j = da(gen);
          if (i == j) continue;
          
          u = A[std::min(i,j)];
          v = A[std::max(i,j)];
          
        } else {
          
          u = A[da(gen)];
          v = B[db(gen)];
        }
        
        uint64_t code = pack_uv(u,v);
        
        // Insert only if not duplicate
        if (seen.insert(code).second)
          local_edges.emplace_back(u,v);
      }
    }
  }
  
  // Merge thread-local edges into global structure
#pragma omp critical
{
  for (auto &e : local_edges) {
    edges.push_back(e);
    
    // Update adjacency
    adj_all[e.first].push_back(e.second);
    adj_all[e.second].push_back(e.first);
  }
}
}

/* -------- global counts -------- */

// Sort adjacency once
sort_adj(adj_all, threads);

// Count wedges & triangles
auto all_wc = count_wedges_triangles_sorted(adj_all, threads);

// Prepare result
List result;

result["all"] = List::create(
  _["edges"] = (long long)edges.size(),
  _["wedges"] = all_wc.first,
  _["triangles"] = all_wc.second
);

/* -------- per W -------- */

// For each probability in probs
for (int pi = 0; pi < probs.size(); ++pi) {
  
  std::mt19937 gen(base_seed + 10007 + pi*997);
  std::bernoulli_distribution bern(probs[pi]);
  
  // Generate W vector
  std::vector<char> W(N);
  for (int i = 0; i < N; ++i)
    W[i] = bern(gen);
  
  // Edge subsets
  std::vector<std::pair<int,int>> edges1, edges2;
  edges1.reserve(edges.size()/4);
  edges2.reserve(edges.size()/2);
  
#pragma omp parallel
{
  std::vector<std::pair<int,int>> le1, le2;
  le1.reserve(1024);
  le2.reserve(1024);
  
#pragma omp for schedule(static)
  for (size_t e = 0; e < edges.size(); ++e) {
    
    int u = edges[e].first;
    int v = edges[e].second;
    
    char wu = W[u];
    char wv = W[v];
    
    // Both endpoints selected
    if ((wu & wv) == 1)
      le1.emplace_back(u,v);
    
    // At least one endpoint selected
    if ((wu | wv) == 1)
      le2.emplace_back(u,v);
  }
  
#pragma omp critical
{
  edges1.insert(edges1.end(), le1.begin(), le1.end());
  edges2.insert(edges2.end(), le2.begin(), le2.end());
}
}

// Helper lambda: build adjacency from edge list
auto build_adj = [&](const std::vector<std::pair<int,int>> &E) {
  
  std::vector<std::vector<int>> adj(N);
  
  for (auto &e : E) {
    adj[e.first].push_back(e.second);
    adj[e.second].push_back(e.first);
  }
  
  return adj;
};

// Build adjacency lists
auto adj1 = build_adj(edges1);
auto adj2 = build_adj(edges2);

// Sort once
sort_adj(adj1, threads);
sort_adj(adj2, threads);

// Count wedges & triangles
auto wc1 = count_wedges_triangles_sorted(adj1, threads);
auto wc2 = count_wedges_triangles_sorted(adj2, threads);

// Store results
std::string pname = "p" + std::to_string(pi+1);

result[pname] = List::create(
  _["edges1"] = (long long)edges1.size(),
  _["wedges1"] = wc1.first,
  _["triangles1"] = wc1.second,
  _["edges2"] = (long long)edges2.size(),
  _["wedges2"] = wc2.first,
  _["triangles2"] = wc2.second
);
}

// Return final result
return result;
}