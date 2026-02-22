// Enable OpenMP support
// [[Rcpp::plugins(openmp)]]

#include <Rcpp.h>          // Rcpp interface
#include <random>          // Random number generators
#include <vector>          // std::vector container
#include <unordered_set>   // Hash set for unique sampling
#include <cinttypes>       // Fixed-width integer types
#include <omp.h>           // OpenMP parallelization
#include <algorithm>       // sort, shuffle
#include <tuple>           // std::pair utilities

using namespace Rcpp;      // Avoid writing Rcpp::

// ------------------------------------------------------------
// Helper: pack (u,v) into one 64-bit integer
// Assumes u,v < 2^32
// ------------------------------------------------------------
static inline uint64_t pack_uv(uint32_t u, uint32_t v) {
  // Place u in high 32 bits and v in low 32 bits
  return ((uint64_t)u << 32) | (uint64_t)v;
}

// ------------------------------------------------------------
// Count wedges and triangles from adjacency list
// ------------------------------------------------------------
std::pair<long long,long long>
  count_wedges_triangles(std::vector<std::vector<int>> &adj, int threads) {
    
    // Number of nodes
    int N = (int)adj.size();
    
    // Initialize counters
    long long wedge_count = 0;
    long long triangle_count = 0;
    
    // Set number of OpenMP threads
    omp_set_num_threads(threads);
    
    // ---- Wedges: sum choose(deg_i, 2) ----
#pragma omp parallel for reduction(+:wedge_count) schedule(static)
    for (int i = 0; i < N; ++i) {
      long long d = (long long)adj[i].size(); // degree
      if (d >= 2)
        wedge_count += (d * (d - 1)) / 2;     // choose(d,2)
    }
    
    // ---- Sort adjacency lists for triangle intersection ----
#pragma omp parallel for schedule(dynamic, 256)
    for (int i = 0; i < N; ++i)
      std::sort(adj[i].begin(), adj[i].end());
    
    // ---- Triangle counting via sorted list intersection ----
#pragma omp parallel for reduction(+:triangle_count) schedule(dynamic, 50)
    for (int i = 0; i < N; ++i) {
      
      const auto &Ni = adj[i]; // neighbors of i
      
      for (size_t idx = 0; idx < Ni.size(); ++idx) {
        
        int j = Ni[idx];
        if (j <= i) continue; // enforce i < j
        
        const auto &Nj = adj[j]; // neighbors of j
        
        size_t a = 0, b = 0;
        
        // Two-pointer merge intersection
        while (a < Ni.size() && b < Nj.size()) {
          
          if (Ni[a] == Nj[b]) {
            int k = Ni[a];
            if (k > j) triangle_count++; // enforce j < k
            a++; b++;
          }
          else if (Ni[a] < Nj[b]) ++a;
          else ++b;
        }
      }
    }
    
    // Return pair (wedges, triangles)
    return {wedge_count, triangle_count};
  }

// ------------------------------------------------------------
// Degree-degree assortativity (Newman Pearson correlation)
// ------------------------------------------------------------
double degree_assortativity_edge_based(
    const std::vector<std::vector<int>> &adj,
    int threads) {
  
  int N = (int)adj.size();
  
  // Compute degrees
  std::vector<int> deg(N);
  for (int i = 0; i < N; ++i)
    deg[i] = (int)adj[i].size();
  
  long long M = 0;              // number of undirected edges
  long long sum_AB = 0;         // sum deg(i)*deg(j)
  long long sum_AplusB = 0;     // sum deg(i)+deg(j)
  long long sum_A2plusB2 = 0;   // sum deg(i)^2 + deg(j)^2
  
  omp_set_num_threads(threads);
  
#pragma omp parallel for reduction(+:M,sum_AB,sum_AplusB,sum_A2plusB2)
  for (int i = 0; i < N; ++i) {
    
    const auto &Ni = adj[i];
    
    for (size_t idx = 0; idx < Ni.size(); ++idx) {
      
      int j = Ni[idx];
      if (j <= i) continue; // count each edge once
      
      int di = deg[i];
      int dj = deg[j];
      
      M += 1;
      sum_AB += (long long)di * dj;
      sum_AplusB += (long long)di + dj;
      sum_A2plusB2 += (long long)di * di +
        (long long)dj * dj;
    }
  }
  
  if (M == 0) return NAN; // no edges
  
  double Md = (double)M;
  
  double E_AB = sum_AB / Md;
  double E_half = (sum_AplusB / Md) * 0.5;
  double E_A2plusB2_over2 =
    (sum_A2plusB2 / Md) * 0.5;
  
  double denom =
    E_A2plusB2_over2 - E_half * E_half;
  
  if (denom == 0.0) return NAN;
  
  double rho =
    (E_AB - E_half * E_half) / denom;
  
  return rho;
}

// ------------------------------------------------------------
// Main: SBM generation + multiple W processing
// ------------------------------------------------------------
// [[Rcpp::export]]
List generate_sbm_edges_with_multiple_W_fast(
    IntegerVector alpha,
    NumericMatrix Pi,
    NumericVector probs,
    int threads = 4,
    int seed = -1) {
  
  const int N = alpha.size();   // number of nodes
  const int K = Pi.nrow();      // number of blocks
  
  if ((int)Pi.ncol() != K)
    Rcpp::stop("Pi must be K x K");
  
  if ((int)probs.size() < 1)
    Rcpp::stop("probs must be length >= 1");
  
  omp_set_num_threads(threads);
  
  // ---- Group nodes by block ----
  std::vector<std::vector<int>> groups(K);
  
  for (int i = 0; i < N; ++i) {
    int b = alpha[i];
    if (b < 1 || b > K)
      Rcpp::stop("alpha values must be in 1..K");
    groups[b-1].push_back(i);
  }
  
  // ---- Base RNG seed ----
  unsigned int base_seed =
    (seed == -1) ?
    std::random_device{}() :
    (unsigned int)seed;
  
  // ---- Store generated edges ----
  std::vector<std::pair<int,int>> edges;
  edges.reserve(1000000);
  
  // ---- Adjacency list for full graph ----
  std::vector<std::vector<int>> adj_all(N);
  
#pragma omp parallel
{
  int tid = omp_get_thread_num();
  std::mt19937 gen(base_seed + 1315423911u + tid*9973);
  
  std::vector<std::pair<int,int>> local_edges;
  
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
      
      uint64_t total_pairs =
        (a == b) ?
        (uint64_t)na * (na - 1) / 2ULL :
        (uint64_t)na * nb;
      
      uint64_t M = 0;
      
      if (total_pairs < (1ULL<<31)) {
        std::binomial_distribution<int>
        bsmall((int)total_pairs, p_ab);
        M = bsmall(gen);
      } else {
        std::poisson_distribution<uint64_t>
        pois(total_pairs * p_ab);
        M = pois(gen);
      }
      
      if (M == 0) continue;
      
      std::unordered_set<uint64_t> seen;
      
      std::uniform_int_distribution<uint32_t>
        uda(0, na - 1);
      std::uniform_int_distribution<uint32_t>
        udb(0, nb - 1);
      
      while (seen.size() < M) {
        
        if (a == b) {
          
          uint32_t ii = uda(gen);
          uint32_t jj = uda(gen);
          if (ii == jj) continue;
          
          uint32_t uidx =
            Ga[std::min(ii,jj)];
          uint32_t vidx =
            Ga[std::max(ii,jj)];
          
          uint64_t code =
            pack_uv(uidx, vidx);
          
          if (seen.insert(code).second)
            local_edges.emplace_back(
              (int)uidx, (int)vidx);
          
        } else {
          
          uint32_t ii = uda(gen);
          uint32_t jj = udb(gen);
          
          uint32_t uidx = Ga[ii];
          uint32_t vidx = Gb[jj];
          
          uint64_t code =
            pack_uv(uidx, vidx);
          
          if (seen.insert(code).second)
            local_edges.emplace_back(
              (int)uidx, (int)vidx);
        }
      }
    }
  }
  
#pragma omp critical
{
  for (auto &e : local_edges) {
    edges.emplace_back(e);
    adj_all[e.first].push_back(e.second);
    adj_all[e.second].push_back(e.first);
  }
}
}

size_t E = edges.size();

// ---- Global counts ----
auto wc_tc_all =
  count_wedges_triangles(adj_all, threads);

double assort_all =
  degree_assortativity_edge_based(adj_all, threads);

List result;

result["all"] = List::create(
  _["edges"] = (long long)E,
  _["wedges"] = wc_tc_all.first,
  _["triangles"] = wc_tc_all.second,
  _["assortativity"] = assort_all
);

// ---- Process multiple W ----
for (int pi = 0; pi < (int)probs.size(); ++pi) {
  
  double prob = probs[pi];
  
  std::mt19937 gen_w(
      base_seed + 10007 + pi*7919);
  
  std::bernoulli_distribution
  bern(prob);
  
  std::vector<char> W(N);
  
  for (int i = 0; i < N; ++i)
    W[i] = bern(gen_w) ? 1 : 0;
  
  std::vector<std::vector<int>> adj1(N), adj2(N);
  
  long long edge_count_1 = 0;
  long long edge_count_2 = 0;
  
  std::vector<int> from_vec,
  to_vec,
  from_block_vec,
  to_block_vec;
  
#pragma omp parallel
{
  std::vector<std::vector<int>> local_adj1(N),
  local_adj2(N);
  
  std::vector<int> local_from,
  local_to,
  local_from_block,
  local_to_block;
  
  long long local_e1 = 0,
    local_e2 = 0;
  
#pragma omp for schedule(static)
  for (size_t ei = 0; ei < E; ++ei) {
    
    int u = edges[ei].first;
    int v = edges[ei].second;
    
    char wu = W[u];
    char wv = W[v];
    
    if ((wu & wv) == 1) {
      
      local_adj1[u].push_back(v);
      local_adj1[v].push_back(u);
      
      local_from.push_back(u+1);
      local_to.push_back(v+1);
      
      local_from_block.push_back(alpha[u]);
      local_to_block.push_back(alpha[v]);
      
      local_e1++;
    }
    
    if ((wu | wv) == 1) {
      
      local_adj2[u].push_back(v);
      local_adj2[v].push_back(u);
      
      local_e2++;
    }
  }
  
#pragma omp critical
{
  for (int i = 0; i < N; ++i) {
    
    adj1[i].insert(adj1[i].end(),
                   local_adj1[i].begin(),
                   local_adj1[i].end());
    
    adj2[i].insert(adj2[i].end(),
                   local_adj2[i].begin(),
                   local_adj2[i].end());
  }
  
  from_vec.insert(from_vec.end(),
                  local_from.begin(),
                  local_from.end());
  
  to_vec.insert(to_vec.end(),
                local_to.begin(),
                local_to.end());
  
  from_block_vec.insert(from_block_vec.end(),
                        local_from_block.begin(),
                        local_from_block.end());
  
  to_block_vec.insert(to_block_vec.end(),
                      local_to_block.begin(),
                      local_to_block.end());
  
  edge_count_1 += local_e1;
  edge_count_2 += local_e2;
}
}

auto wc_tc_1 =
  count_wedges_triangles(adj1, threads);

auto wc_tc_2 =
  count_wedges_triangles(adj2, threads);

std::string pname =
  "p" + std::to_string(pi+1);

result[pname] = List::create(
  _["edges1"] = edge_count_1,
  _["wedges1"] = wc_tc_1.first,
  _["triangles1"] = wc_tc_1.second,
  _["edge_list1"] = DataFrame::create(
    _["from"] = from_vec,
    _["to"] = to_vec,
    _["from_block"] = from_block_vec,
    _["to_block"] = to_block_vec
  ),
  _["edges2"] = edge_count_2,
  _["wedges2"] = wc_tc_2.first,
  _["triangles2"] = wc_tc_2.second
);
}

return result;
}
