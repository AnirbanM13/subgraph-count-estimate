library(Rcpp)
library(parallel)

# =========================
# generate_Rt_cpp
# =========================
cppFunction('
#include <Rcpp.h>
using namespace Rcpp;

void combn_recursive(int offset, int k, int n,
                     std::vector<int>& current,
                     std::vector< std::vector<int> >& result) {
  if (current.size() == k) {
    result.push_back(current);
    return;
  }
  for (int i = offset; i <= n; i++) {
    current.push_back(i);
    combn_recursive(i + 1, k, n, current, result);
    current.pop_back();
  }
}

std::vector< std::vector<int> > generate_permutations(int t) {
  std::vector<int> perm(t);
  for (int i = 0; i < t; i++) perm[i] = i + 1;
  
  std::vector< std::vector<int> > perms;
  do {
    perms.push_back(perm);
  } while (std::next_permutation(perm.begin(), perm.end()));
  
  return perms;
}

// [[Rcpp::export]]
List generate_Rt_cpp(int R, int t) {
  
  std::vector< std::vector<int> > combs;
  std::vector<int> current;
  combn_recursive(1, t, R, current, combs);
  
  std::vector< std::vector<int> > perms = generate_permutations(t);
  
  int total = combs.size() * combs.size() * perms.size();
  List Rt(total);
  
  int idx = 0;
  
  for (size_t i = 0; i < combs.size(); i++) {
    for (size_t j = 0; j < combs.size(); j++) {
      for (size_t p = 0; p < perms.size(); p++) {
        
        Rt[idx] = List::create(
          Named("j") = combs[i],
          Named("k") = combs[j],
          Named("xi") = perms[p]
        );
        
        idx++;
      }
    }
  }
  
  return Rt;
}
')


# =========================
# theta_model_sum_cpp (B term)
# =========================
cppFunction('
#include <Rcpp.h>
using namespace Rcpp;

void dfs(int pos,
         int R, int K, int t,
         IntegerVector j,
         IntegerVector k,
         IntegerVector xi,
         NumericMatrix Pi,
         NumericVector lambda,
         std::vector<int>& u,
         std::vector<int>& v,
         std::vector< std::pair<int,int> >& edges_H,
         std::vector< std::pair<int,int> >& edges_H_minus_B,
         double& total) {
  
  if (pos == 2*R) {
    
    double prod_u = 1.0;
    for (auto &e : edges_H)
      prod_u *= Pi(u[e.first], u[e.second]);
    
    double prod_v = 1.0;
    for (auto &e : edges_H_minus_B)
      prod_v *= Pi(v[e.first], v[e.second]);
    
    double prod_lambda_u = 1.0;
    for (int i = 0; i < R; i++)
      prod_lambda_u *= lambda[u[i]];
    
    double prod_lambda_v = 1.0;
    for (int i = 0; i < R; i++) {
      bool inAk = false;
      for (int r = 0; r < t; r++) {
        if (k[r]-1 == i) {
          inAk = true;
          break;
        }
      }
      if (!inAk)
        prod_lambda_v *= lambda[v[i]];
    }
    
    total += prod_u * prod_v * prod_lambda_u * prod_lambda_v;
    return;
  }
  
  if (pos < R) {
    for (int val = 0; val < K; val++) {
      u[pos] = val;
      dfs(pos+1, R, K, t, j, k, xi, Pi, lambda,
          u, v, edges_H, edges_H_minus_B, total);
    }
  } else {
    int idx = pos - R;
    
    bool constrained = false;
    int forced_val = -1;
    
    for (int r = 0; r < t; r++) {
      if (k[xi[r]-1] - 1 == idx) {
        constrained = true;
        forced_val = u[j[r]-1];
        break;
      }
    }
    
    if (constrained) {
      v[idx] = forced_val;
      dfs(pos+1, R, K, t, j, k, xi, Pi, lambda,
          u, v, edges_H, edges_H_minus_B, total);
    } else {
      for (int val = 0; val < K; val++) {
        v[idx] = val;
        dfs(pos+1, R, K, t, j, k, xi, Pi, lambda,
            u, v, edges_H, edges_H_minus_B, total);
      }
    }
  }
}


// [[Rcpp::export]]
double theta_model_sum_cpp(IntegerVector j,
                             IntegerVector k,
                             IntegerVector xi,
                             NumericMatrix Pi,
                             NumericMatrix adj,
                             NumericVector lambda) {
  
  int R = adj.nrow();
  int K = lambda.size();
  int t = j.size();
  
  std::vector< std::pair<int,int> > edges_H;
  for (int i = 0; i < R; i++) {
    for (int j2 = i+1; j2 < R; j2++) {
      if (adj(i,j2) == 1)
        edges_H.push_back({i,j2});
    }
  }
  
  std::vector<bool> inAk(R, false);
  for (int i = 0; i < t; i++)
    inAk[k[i]-1] = true;
  
  std::vector< std::pair<int,int> > edges_H_minus_B;
  for (auto &e : edges_H) {
    if (!(inAk[e.first] && inAk[e.second]))
      edges_H_minus_B.push_back(e);
  }
  
  std::vector<int> u(R, 0), v(R, 0);
  double total = 0.0;
  
  dfs(0, R, K, t, j, k, xi, Pi, lambda,
      u, v, edges_H, edges_H_minus_B, total);
  
  return total;
}
')


# =========================
# vertex cover pieces (A term)
# =========================
H1_union_H2 <- function(R, t, edges_H, j, k, xi) {
  
  parent <- 1:(2*R)
  
  find <- function(x) {
    if (parent[x] != x) parent[x] <<- Recall(parent[x])
    parent[x]
  }
  
  unite <- function(x, y) {
    px <- find(x); py <- find(y)
    if (px != py) parent[px] <<- py
  }
  
  for (r in seq_len(t)) {
    unite(j[r], k[xi[r]] + R)
  }
  
  edge_set <- list()
  add_edge <- function(u, v) {
    if (u == v) return()
    e <- sort(c(u, v))
    edge_set[[paste(e, collapse = "-")]] <<- e
  }
  
  for (e in edges_H) {
    u <- e[1]; v <- e[2]
    add_edge(find(u), find(v))
    add_edge(find(u + R), find(v + R))
  }
  
  edges <- do.call(rbind, edge_set)
  if (is.null(edges)) edges <- matrix(numeric(0), ncol = 2)
  
  verts <- unique(as.vector(edges))
  id_map <- setNames(seq_along(verts), verts)
  
  edges_comp <- t(apply(edges, 1, function(e) {
    c(id_map[as.character(e[1])],
      id_map[as.character(e[2])])
  }))
  
  V <- length(verts)
  
  is_vertex_cover <- function(mask) {
    for (i in seq_len(nrow(edges_comp))) {
      u <- edges_comp[i,1]; v <- edges_comp[i,2]
      if (!bitwAnd(mask, bitwShiftL(1, u-1)) &&
          !bitwAnd(mask, bitwShiftL(1, v-1)))
        return(FALSE)
    }
    TRUE
  }
  
  all_vertex_covers <- list()
  idx <- 1
  
  for (mask in 0:(2^V - 1)) {
    if (is_vertex_cover(mask)) {
      bits <- as.integer(intToBits(mask))[1:V]
      all_vertex_covers[[idx]] <- which(bits == 1)
      idx <- idx + 1
    }
  }
  
  list(all_vertex_covers = all_vertex_covers)
}


psi_design_cov <- function(covers, R, t, p) {
  sizes <- sapply(covers, length)
  sum(p^sizes * (1 - p)^(2*R - t - sizes))
}


get_all_vertex_covers <- function(adj_matrix) {
  R <- nrow(adj_matrix)
  
  edges_idx <- which(adj_matrix == 1, arr.ind = TRUE)
  edges_idx <- edges_idx[edges_idx[,1] < edges_idx[,2], , drop = FALSE]
  edges <- split(edges_idx, seq(nrow(edges_idx)))
  edges <- lapply(edges, function(e) c(e[1], e[2]))
  
  covers <- list(); idx <- 1
  
  for (mask in 0:(2^R - 1)) {
    bits <- as.integer(intToBits(mask))[1:R]
    ok <- TRUE
    for (e in edges) {
      if (bits[e[1]] == 0 && bits[e[2]] == 0) {
        ok <- FALSE; break
      }
    }
    if (ok) {
      covers[[idx]] <- which(bits == 1)
      idx <- idx + 1
    }
  }
  
  covers
}


f2_vertex_cover <- function(adj_matrix, p) {
  covers <- get_all_vertex_covers(adj_matrix)
  sizes <- sapply(covers, length)
  sum(p^sizes * (1 - p)^(nrow(adj_matrix) - sizes))
}


# =========================
# FINAL: var_H
# =========================
var_H <- function(
    adj_matrix, Pi, lambda,
    N, pN,
    mode = c("ego-centric", "induced"),
    mc.cores = 8
) {
  
  mode <- match.arg(mode)
  R <- nrow(adj_matrix)
  
  edges_H <- which(adj_matrix == 1, arr.ind = TRUE)
  edges_H <- edges_H[edges_H[,1] < edges_H[,2], , drop = FALSE]
  edges_H <- split(edges_H, seq(nrow(edges_H)))
  edges_H <- lapply(edges_H, function(e) c(e[1], e[2]))
  
  total_sum <- 0
  
  for (t in 1:R) {
    
    Rt_list <- generate_Rt_cpp(R, t)
    
    vals <- mclapply(Rt_list, function(triple) {
      
      j  <- triple$j
      k  <- triple$k
      xi <- triple$xi
      
      B_val <- theta_model_sum_cpp(
        j, k, xi, Pi, adj_matrix, lambda
      )
      
      if (mode == "induced") return(B_val)
      
      res <- H1_union_H2(R, t, edges_H, j, k, xi)
      
      A_val <- psi_design_cov(
        res$all_vertex_covers, R, t, pN
      )
      
      (A_val - (f2_vertex_cover(adj_matrix, pN))^2) * B_val
      
    }, mc.cores = mc.cores)
    
    sum_t <- sum(unlist(vals))
    
    if (mode == "induced") {
      total_sum <- total_sum +
        sum_t * (N^(2*R - t)) * (pN^(2*R - t)*(1-pN^t))
    } else {
      total_sum <- total_sum +
        sum_t * (N^(2*R - t))
    }
  }
  
  total_sum
}




# tri <- matrix(
#   c(0,1,1,
#     1,0,1,
#     1,1,0), nrow = 3, byrow = T)
# 
# 
# 
# 
# N.main <- 20000
# al.vec <- c(rep(1,N.main/2), rep(2,N.main/2))
# Pi.mat <- rbind(c(0.2, 0.05), c(0.05, 0.1))
# 
# p.vec = c(1,5)/100
# #M = 2000
# lambda <- c(0.5,0.5)
# 
# 
# 
# 
# 
# var_new_tri <- var_H(
#   tri, Pi.mat, lambda,
#   N = N.main,
#   pN = 0.01,
#   mode = "ego-centric",
#   mc.cores = 4
# )
# 
# var_new_tri/ var_Triangle_ego(N.main, al.vec, Pi.mat, 0.01)















