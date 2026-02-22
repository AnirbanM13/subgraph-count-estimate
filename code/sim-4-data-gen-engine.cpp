// Enable OpenMP parallelization support
// [[Rcpp::plugins(openmp)]]

#include <Rcpp.h>        // Rcpp interface
#include <random>        // Random number generators
#include <vector>        // std::vector container
#include <cinttypes>     // Fixed-width integer types
#include <omp.h>         // OpenMP parallelization
#include <algorithm>     // std::swap, std::round
#include <unordered_set> // Hash set for uniqueness control

using namespace Rcpp;    // Avoid writing Rcpp::


/* ============================================================
 PCG32 RANDOM NUMBER GENERATOR
 ============================================================ */

// Define PCG32 RNG structure
struct pcg32 {
  
  using result_type = uint32_t;   // Required typedef for RNG compatibility
  
  uint64_t state;  // Internal state
  uint64_t inc;    // Stream increment
  
  // Constructor with seed and sequence
  pcg32(uint64_t seed = 42u, uint64_t seq = 54u) {
    seed_rng(seed, seq);
  }
  
  // Seed initialization
  void seed_rng(uint64_t seed, uint64_t seq) {
    state = 0u;
    inc = (seq << 1u) | 1u;  // Make increment odd
    operator()();            // Advance state
    state += seed;
    operator()();            // Advance again
  }
  
  // RNG step
  inline uint32_t operator()() {
    uint64_t oldstate = state;
    
    // Linear congruential update
    state = oldstate * 6364136223846793005ULL + inc;
    
    // Output transformation
    uint32_t xorshifted =
      ((oldstate >> 18u) ^ oldstate) >> 27u;
    
    uint32_t rot = oldstate >> 59u;
    
    // Rotate bits
    return (xorshifted >> rot) |
      (xorshifted << ((-rot) & 31));
  }
  
  // Required by <random>
  static constexpr uint32_t min() { return 0u; }
  static constexpr uint32_t max() { return UINT32_MAX; }
};


/* ============================================================
 Helper: pack (u,v) into 64-bit integer
 ============================================================ */

// Combines two 32-bit integers into one 64-bit value
static inline uint64_t pack_uv(uint32_t u, uint32_t v) {
  return ((uint64_t)u << 32) | (uint64_t)v;
}


/* ============================================================
 MAIN FUNCTION
 ============================================================ */

// Export function to R
// [[Rcpp::export]]
List generate_sbm_edges_streamed_fast(
    IntegerVector alpha,    // Block labels
    NumericMatrix Pi,       // Block probability matrix
    NumericVector probs,    // W probabilities
    int threads = 4,        // Number of OpenMP threads
    int seed = -1) {        // Random seed
  
  // Number of nodes
  const int N = alpha.size();
  
  // Number of blocks
  const int K = Pi.nrow();
  
  // Sanity checks
  if ((int)Pi.ncol() != K)
    Rcpp::stop("Pi must be K x K");
  
  if ((int)probs.size() < 1)
    Rcpp::stop("probs must be length >= 1");
  
  // Set number of threads
  omp_set_num_threads(threads);
  
  /* ------------------------
   Group nodes by block
   ------------------------ */
  
  std::vector<std::vector<int>> groups(K);
  
  for (int i = 0; i < N; ++i) {
    int b = alpha[i];
    if (b < 1 || b > K)
      Rcpp::stop("alpha values must be in 1..K");
    
    groups[b-1].push_back(i);
  }
  
  // Base seed
  unsigned int base_seed =
    (seed == -1) ?
    std::random_device{}() :
    (unsigned int)seed;
  
  /* ------------------------
   Thread-local edge storage
   ------------------------ */
  
  std::vector<std::vector<uint32_t>> edges_u_local(threads);
  std::vector<std::vector<uint32_t>> edges_v_local(threads);
  
  
#pragma omp parallel
{
  int tid = omp_get_thread_num();
  
  // Thread-specific RNG
  pcg32 gen(base_seed + 1315423911u + tid*9973, tid);
  
  // References to thread-local buffers
  auto &eu = edges_u_local[tid];
  auto &ev = edges_v_local[tid];
  
  eu.reserve(100000);
  ev.reserve(100000);
  
#pragma omp for schedule(dynamic)
  for (int a = 0; a < K; ++a) {
    for (int b = a; b < K; ++b) {
      
      const auto &Ga = groups[a];
      const auto &Gb = groups[b];
      
      int na = Ga.size();
      int nb = Gb.size();
      
      if (na == 0 || nb == 0) continue;
      
      double p_ab = Pi(a,b);
      if (p_ab <= 0.0) continue;
      
      // Total possible pairs
      uint64_t total_pairs =
        (a == b)
        ? (uint64_t)na * (na - 1) / 2ULL
      : (uint64_t)na * nb;
      
      uint64_t M = 0;
      
      // Binomial sampling
      if (total_pairs < (1ULL<<31)) {
        std::binomial_distribution<int>
        bsmall((int)total_pairs, p_ab);
        M = bsmall(gen);
      } else {
        if (p_ab < 1e-3) {
          std::poisson_distribution<uint64_t>
          pois(total_pairs * p_ab);
          M = pois(gen);
        } else {
          M = std::round(total_pairs * p_ab);
          if (M > total_pairs) M = total_pairs;
        }
      }
      
      if (M == 0) continue;
      
      /* ------------------------
       Sampling strategies
       ------------------------ */
      
      // Sparse regime
      if (M * 100 < total_pairs) {
        
        std::uniform_int_distribution<uint32_t>
        uda(0, na - 1);
        std::uniform_int_distribution<uint32_t>
          udb(0, nb - 1);
        
        for (uint64_t k = 0; k < M; ++k) {
          
          uint32_t uidx, vidx;
          
          if (a == b) {
            do {
              uidx = Ga[uda(gen)];
              vidx = Ga[uda(gen)];
            } while (uidx == vidx);
            
            if (uidx > vidx)
              std::swap(uidx, vidx);
            
          } else {
            uidx = Ga[uda(gen)];
            vidx = Gb[udb(gen)];
          }
          
          eu.push_back(uidx);
          ev.push_back(vidx);
        }
        
        // Dense regime
      } else if (M > total_pairs * 0.9) {
        
        uint64_t missing = total_pairs - M;
        
        std::unordered_set<uint64_t> skip;
        skip.reserve(missing*2 + 10);
        
        std::uniform_int_distribution<uint32_t>
          uda(0, na - 1);
        std::uniform_int_distribution<uint32_t>
          udb(0, nb - 1);
        
        while (skip.size() < missing) {
          
          uint32_t uidx, vidx;
          
          if (a == b) {
            do {
              uidx = Ga[uda(gen)];
              vidx = Ga[uda(gen)];
            } while (uidx == vidx);
            
            if (uidx > vidx)
              std::swap(uidx, vidx);
          } else {
            uidx = Ga[uda(gen)];
            vidx = Gb[udb(gen)];
          }
          
          skip.insert(pack_uv(uidx, vidx));
        }
        
        for (int ii : Ga) {
          for (int jj : Gb) {
            
            if (a == b && ii >= jj) continue;
            
            if (skip.count(pack_uv(ii,jj)) == 0) {
              eu.push_back(ii);
              ev.push_back(jj);
            }
          }
        }
        
        // Moderate density regime
      } else {
        
        std::unordered_set<uint64_t> seen;
        seen.reserve(M*2 + 10);
        
        std::uniform_int_distribution<uint32_t>
          uda(0, na - 1);
        std::uniform_int_distribution<uint32_t>
          udb(0, nb - 1);
        
        while (seen.size() < M) {
          
          uint32_t uidx, vidx;
          
          if (a == b) {
            do {
              uidx = Ga[uda(gen)];
              vidx = Ga[uda(gen)];
            } while (uidx == vidx);
            
            if (uidx > vidx)
              std::swap(uidx, vidx);
          } else {
            uidx = Ga[uda(gen)];
            vidx = Gb[udb(gen)];
          }
          
          uint64_t code =
            pack_uv(uidx, vidx);
          
          if (seen.insert(code).second) {
            eu.push_back(uidx);
            ev.push_back(vidx);
          }
        }
      }
    }
  }
} // End OpenMP region


/* ------------------------
 Merge thread-local edges
 ------------------------ */

std::vector<uint32_t> edges_u, edges_v;

size_t total_E = 0;

for (int t = 0; t < threads; ++t)
  total_E += edges_u_local[t].size();

edges_u.reserve(total_E);
edges_v.reserve(total_E);

for (int t = 0; t < threads; ++t) {
  edges_u.insert(edges_u.end(),
                 edges_u_local[t].begin(),
                 edges_u_local[t].end());
  
  edges_v.insert(edges_v.end(),
                 edges_v_local[t].begin(),
                 edges_v_local[t].end());
}

size_t E = edges_u.size();


/* ------------------------
 Precompute W for all probs
 ------------------------ */

std::vector<std::vector<char>> Wmat(
    probs.size(),
    std::vector<char>(N));

for (int pi = 0; pi < probs.size(); ++pi) {
  
  pcg32 gen_w(base_seed + 10007 + pi*7919, pi);
  
  std::bernoulli_distribution
  bern(probs[pi]);
  
  for (int i = 0; i < N; ++i)
    Wmat[pi][i] =
      bern(gen_w) ? 1 : 0;
}


/* ------------------------
 Count edges under W
 ------------------------ */

List result;

result["all"] =
  List::create(
    _["edges"] = (long long)E);

for (int pi = 0; pi < probs.size(); ++pi) {
  
  long long edge_count_1 = 0;
  long long edge_count_2 = 0;
  
  const auto &W = Wmat[pi];
  
#pragma omp parallel for reduction(+:edge_count_1,edge_count_2)
  for (size_t ei = 0; ei < E; ++ei) {
    
    char wu = W[edges_u[ei]];
    char wv = W[edges_v[ei]];
    
    if ((wu & wv) == 1)
      edge_count_1++;
    
    if ((wu | wv) == 1)
      edge_count_2++;
  }
  
  std::string pname =
    "p" + std::to_string(pi+1);
  
  result[pname] =
    List::create(
      _["edges1"] = edge_count_1,
      _["edges2"] = edge_count_2
    );
}

return result;
}