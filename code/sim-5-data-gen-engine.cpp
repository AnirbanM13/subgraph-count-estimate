// Enable OpenMP support in Rcpp (for parallel execution)
// [[Rcpp::plugins(openmp)]]

#include <Rcpp.h>      // Rcpp interface between R and C++
#include <random>      // C++ random number distributions
#include <vector>      // std::vector container
#include <omp.h>       // OpenMP header for multithreading

using namespace Rcpp;  // Avoid writing Rcpp:: everywhere

// =====================================================
//                PCG32 RANDOM NUMBER GENERATOR
// =====================================================
// This is a fast, statistically strong 32-bit RNG.
// It is standard-compliant with <random> requirements.
struct pcg32 {
  
  using result_type = uint32_t;   // Required typedef for C++ RNG engines
  
  uint64_t state;  // Internal RNG state
  uint64_t inc;    // Stream selector (must be odd)
  
  // Constructor with optional seed and sequence
  pcg32(uint64_t seed = 42u, uint64_t seq = 54u) {
    seed_rng(seed, seq);  // Initialize RNG state
  }
  
  // Seed initialization function
  void seed_rng(uint64_t seed, uint64_t seq) {
    state = 0u;                     // Reset state
    inc = (seq << 1u) | 1u;         // Make sure increment is odd
    operator()();                   // Advance state
    state += seed;                  // Add seed
    operator()();                   // Advance again
  }
  
  // Core RNG function
  inline uint32_t operator()() {
    uint64_t oldstate = state;  // Save current state
    
    // Update internal state (LCG step)
    state = oldstate * 6364136223846793005ULL + inc;
    
    // Output transformation (xorshift + rotate)
    uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
    uint32_t rot = oldstate >> 59u;
    
    // Rotate bits to produce final value
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
  }
  
  // Required by C++ random engine concept
  result_type min() const noexcept { return 0u; }
  result_type max() const noexcept { return UINT32_MAX; }
};

// =====================================================
//         MAIN FUNCTION EXPORTED TO R
// =====================================================
// [[Rcpp::export]]
List generate_sbm_edges_fast_no_edgelist(
    IntegerVector alpha,   // Node block memberships (size N)
    NumericMatrix Pi,      // K x K SBM probability matrix
    NumericVector probs,   // Vector of Bernoulli probabilities
    int threads = 4,       // Number of OpenMP threads
    int seed = -1) {       // Optional seed
  
  const int N = alpha.size();   // Number of nodes
  const int K = Pi.nrow();      // Number of blocks
  const int P = probs.size();   // Number of probability settings
  
  // Ensure Pi is square
  if (Pi.ncol() != K)
    stop("Pi must be K x K");
  
  // Set number of threads for OpenMP
  omp_set_num_threads(threads);
  
  // -------------------------------------------------
  // Compute number of nodes in each block
  // -------------------------------------------------
  std::vector<int64_t> n_block(K, 0);  // Initialize block counts to 0
  
  for (int i = 0; i < N; ++i) {
    int b = alpha[i] - 1;   // Convert R 1-based index to C++ 0-based
    
    if (b < 0 || b >= K)
      stop("alpha values must be in 1..K");
    
    n_block[b]++;  // Increment block count
  }
  
  // -------------------------------------------------
  // Compute W counts per block (for each probability)
  // -------------------------------------------------
  std::vector<std::vector<int64_t>> n1(P, std::vector<int64_t>(K, 0));
  std::vector<std::vector<int64_t>> n0(P, std::vector<int64_t>(K, 0));
  // n1[p][b] = number of nodes in block b with W=1
  // n0[p][b] = number of nodes in block b with W=0
  
  // Base seed (random if seed == -1)
  unsigned int base_seed =
    (seed == -1) ? std::random_device{}() : (unsigned int)seed;
  
  // Parallel over probabilities
#pragma omp parallel for
  for (int p = 0; p < P; ++p) {
    
    // Independent RNG per probability
    pcg32 gen(base_seed + 10007 + p * 7919, p);
    
    // Bernoulli distribution with probability probs[p]
    std::bernoulli_distribution bern(probs[p]);
    
    // For each node
    for (int i = 0; i < N; ++i) {
      int b = alpha[i] - 1;   // Block index
      
      if (bern(gen))
        n1[p][b]++;  // Node assigned W=1
      else
        n0[p][b]++;  // Node assigned W=0
    }
  }
  
  // -------------------------------------------------
  // Compute total number of SBM edges (no W)
  // -------------------------------------------------
  int64_t total_edges = 0;
  
  {
    pcg32 gen(base_seed + 99991, 123);  // RNG for total edges
    
    for (int a = 0; a < K; ++a) {
      for (int b = a; b < K; ++b) {
        
        double pij = Pi(a, b);  // Edge probability
        if (pij <= 0.0) continue;
        
        int64_t na = n_block[a];
        int64_t nb = n_block[b];
        
        // Total possible edges between blocks
        int64_t T = (a == b)
          ? na * (na - 1) / 2   // Within-block
        : na * nb;            // Between-block
        
        if (T <= 0) continue;
        
        // If T small enough, use binomial
        if (T < (1LL << 31)) {
          std::binomial_distribution<int> bin(T, pij);
          total_edges += bin(gen);
        } else {
          // If very large, approximate with Poisson
          std::poisson_distribution<int64_t> pois(T * pij);
          total_edges += pois(gen);
        }
      }
    }
  }
  
  // Store total result
  List result;
  result["all"] = List::create(_["edges"] = total_edges);
  
  // -------------------------------------------------
  // Compute per-probability edge counts
  // -------------------------------------------------
  for (int p = 0; p < P; ++p) {
    
    pcg32 gen(base_seed + 50021 + p * 1777, p + 7);
    
    int64_t edges1 = 0;  // Edges where both W=1
    int64_t edges2 = 0;  // Edges where at least one W=1
    
    for (int a = 0; a < K; ++a) {
      for (int b = a; b < K; ++b) {
        
        double pij = Pi(a, b);
        if (pij <= 0.0) continue;
        
        int64_t na = n_block[a];
        int64_t nb = n_block[b];
        
        int64_t T = (a == b)
          ? na * (na - 1) / 2
        : na * nb;
        
        if (T <= 0) continue;
        
        int64_t M;
        
        // Sample total edges for this block pair
        if (T < (1LL << 31)) {
          std::binomial_distribution<int> bin(T, pij);
          M = bin(gen);
        } else {
          std::poisson_distribution<int64_t> pois(T * pij);
          M = pois(gen);
        }
        
        if (M == 0) continue;
        
        // ---- Edges where both nodes have W=1 ----
        int64_t T11 = (a == b)
          ? n1[p][a] * (n1[p][a] - 1) / 2
        : n1[p][a] * n1[p][b];
        
        if (T11 > 0) {
          std::binomial_distribution<int64_t> bin11(T11, pij);
          edges1 += bin11(gen);
        }
        
        // ---- Edges where both nodes have W=0 ----
        int64_t T00 = (a == b)
          ? n0[p][a] * (n0[p][a] - 1) / 2
        : n0[p][a] * n0[p][b];
        
        int64_t e00 = 0;
        if (T00 > 0) {
          std::binomial_distribution<int64_t> bin00(T00, pij);
          e00 = bin00(gen);
        }
        
        // Edges where at least one W=1
        edges2 += (M - e00);
      }
    }
    
    // Store results per probability
    std::string pname = "p" + std::to_string(p + 1);
    
    result[pname] = List::create(
      _["edges1"] = edges1,
      _["edges2"] = edges2
    );
  }
  
  return result;  // Return to R
}
