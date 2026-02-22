############################################################
# Clean Environment
############################################################
rm(list = ls())

############################################################
# Load Political Blog Network (Edge List)
############################################################

url <- "https://raw.githubusercontent.com/NiccoloBalestrieri/Political-Blog-2004-U.S.-Election-Analysis/main/dataset/edge_list.csv"

edges_raw <- read.csv(url, sep = ";", stringsAsFactors = FALSE)
colnames(edges_raw) <- c("from", "to")

edges_clean <- edges_raw[
  !is.na(edges_raw$from) &
    !is.na(edges_raw$to)   &
    edges_raw$from != ""   &
    edges_raw$to   != "",
]

edges_clean$from <- as.integer(edges_clean$from)
edges_clean$to   <- as.integer(edges_clean$to)

edges_clean <- edges_clean[
  !is.na(edges_clean$from) &
    !is.na(edges_clean$to),
]

edges_clean <- edges_clean[edges_clean$from != edges_clean$to, ]
edges_clean <- unique(edges_clean)

############################################################
# Compile C++ + Load R Files
############################################################

Rcpp::sourceCpp("Desktop/r-file/edge-list-to-adj-mat-conv-1.cpp")
Rcpp::sourceCpp("Desktop/r-file/pi_hat.cpp")

source("~/Desktop/r-file/bethe_hessian_clustering.R")
source("~/Desktop/r-file/var-exp-final.R")

############################################################
# Build Adjacency Matrix
############################################################

A <- edge_list_to_adj_unique(edges_clean, directed = FALSE) # make the network undirected

N <- length(unique(c(edges_clean$from, edges_clean$to)))

############################################################
# Load Orientation Labels
############################################################

url_orient <- "https://raw.githubusercontent.com/NiccoloBalestrieri/Political-Blog-2004-U.S.-Election-Analysis/main/dataset/orientation.csv"

orientation_raw <- scan(url_orient, what = character(), sep = ",", quiet = TRUE)
al.vec <- as.integer(factor(orientation_raw))

############################################################
# Graph Statistics Functions
############################################################

edge_count <- function(A) {
  sum(A)
}

wedge_count <- function(A) {
  deg <- rowSums(A)
  sum(deg * (deg - 1))
}

library(Rcpp)
library(Matrix)

Rcpp::cppFunction('
double triangle_count_rcpp(IntegerVector p, IntegerVector i, int n) {
  double tri = 0.0;

  for (int col = 0; col < n; col++) {
    for (int idx = p[col]; idx < p[col + 1]; idx++) {

      int row = i[idx];
      if (row <= col) continue;

      int a = p[col], b = p[row];
      int a_end = p[col + 1], b_end = p[row + 1];

      while (a < a_end && b < b_end) {
        if (i[a] == i[b]) {
          tri++;
          a++; b++;
        } else if (i[a] < i[b]) {
          a++;
        } else {
          b++;
        }
      }
    }
  }
  return tri;
}
')

triangle_count_fast <- function(A) {
  A <- as(A, "dgCMatrix")
  triangle_count_rcpp(A@p, A@i, ncol(A))
}

CC_cnt <- function(A) {
  W <- wedge_count(A)
  Tri <- triangle_count_fast(A)
  ifelse(W != 0, 2 * Tri / W, 0)
}

############################################################
# Population Statistics
############################################################

E_true  <- edge_count(A)
S_true  <- wedge_count(A)
T_true  <- triangle_count_fast(A)
CC_true <- CC_cnt(A)

############################################################
# Estimate Block Probability Matrix
############################################################

Pi.mat <- PI_hat_cpp(al.vec, A)

############################################################
# Sampling Function
############################################################

sample_induced_egocentric_R <- function(A, p = 0.5) {
  
  if (!is.matrix(A)) stop("A must be a matrix")
  if (nrow(A) != ncol(A)) stop("Adjacency matrix must be square")
  
  N <- nrow(A)
  
  W <- rbinom(N, 1, p)
  S <- which(W == 1)
  
  A_induced <- A[S, S, drop = FALSE]
  
  M_ego <- outer(W, W, pmax)
  A_ego <- A * M_ego
  
  list(
    W = W,
    sampled_nodes = S,
    A_induced = A_induced,
    A_egocentric = A_ego
  )
}

############################################################
# Run Sampling
############################################################

set.seed(8)
p <- 0.1

A_est <- sample_induced_egocentric_R(A, p)

A.ind <- A_est$A_induced
A.ego <- A_est$A_egocentric

com.ind    <- bethe_hessian_clustering_adj(A.ind)
Pi.mat.ind <- PI_hat_cpp(com.ind, A.ind)

############################################################
# Evaluation Function
############################################################

evaluate_sampling <- function(A, A.ind, A.ego,
                              N, com.ind, Pi.mat.ind, p) {
  
  # True values
  E_true  <- edge_count(A)
  S_true  <- wedge_count(A)
  T_true  <- triangle_count_fast(A)
  CC_true <- CC_cnt(A)
  
  # IND
  v_e_i <- var_Edge_ind(N, com.ind, Pi.mat.ind, p)
  v_s_i <- var_wedge_ind(N, com.ind, Pi.mat.ind, p)
  v_t_i <- var_Triangle_ind(N, com.ind, Pi.mat.ind, p)
  v_c_i <- var_cc_ind(N, com.ind, Pi.mat.ind, p, 0)
  
  e_hat <- edge_count(A.ind)
  s_hat <- wedge_count(A.ind)
  t_hat <- triangle_count_fast(A.ind)
  c_hat <- CC_cnt(A.ind)
  
  E_ind_est <- e_hat / p^2
  S_ind_est <- s_hat / p^3
  T_ind_est <- t_hat / p^3
  C_ind_est <- c_hat
  
  lo_e_ind <- (1/p^2) * (e_hat - 1.96 * sqrt(v_e_i))
  up_e_ind <- (1/p^2) * (e_hat + 1.96 * sqrt(v_e_i))
  
  lo_s_ind <- (1/p^3) * (s_hat - 1.96 * sqrt(v_s_i))
  up_s_ind <- (1/p^3) * (s_hat + 1.96 * sqrt(v_s_i))
  
  lo_t_ind <- (1/p^3) * (t_hat - 1.96 * sqrt(v_t_i))
  up_t_ind <- (1/p^3) * (t_hat + 1.96 * sqrt(v_t_i))
  
  lo_c_ind <- c_hat - 1.96 * sqrt(v_c_i)
  up_c_ind <- c_hat + 1.96 * sqrt(v_c_i)
  
  # EGO
  v_e_e <- var_Edge_ego(N, com.ind, Pi.mat.ind, p)
  v_s_e <- var_wedge_ego(N, com.ind, Pi.mat.ind, p)
  v_t_e <- var_Triangle_ego(N, com.ind, Pi.mat.ind, p)
  v_c_e <- var_cc_ego(N, com.ind, Pi.mat.ind, p, 0)
  
  e_hat1 <- edge_count(A.ego)
  s_hat1 <- wedge_count(A.ego)
  t_hat1 <- triangle_count_fast(A.ego)
  c_hat1 <- CC_cnt(A.ego)
  
  E_scale <- 1/(2*p - p^2)
  S_scale <- 1/(p + (1-p)*p^2)
  T_scale <- 1/(3*(1-p)*p^2 + p^3)
  f <- (p + (1-p)*p^2) / (3*p^2*(1-p) + p^3)
  
  E_ego_est <- e_hat1 * E_scale
  S_ego_est <- s_hat1 * S_scale
  T_ego_est <- t_hat1 * T_scale
  C_ego_est <- c_hat1 * f
  
  lo_e_ego <- E_scale * (e_hat1 - 1.96 * sqrt(v_e_e))
  up_e_ego <- E_scale * (e_hat1 + 1.96 * sqrt(v_e_e))
  
  lo_s_ego <- S_scale * (s_hat1 - 1.96 * sqrt(v_s_e))
  up_s_ego <- S_scale * (s_hat1 + 1.96 * sqrt(v_s_e))
  
  lo_t_ego <- T_scale * (t_hat1 - 1.96 * sqrt(v_t_e))
  up_t_ego <- T_scale * (t_hat1 + 1.96 * sqrt(v_t_e))
  
  lo_c_ego <- C_ego_est - 1.96 * sqrt(v_c_e)
  up_c_ego <- C_ego_est + 1.96 * sqrt(v_c_e)
  
  coverage_vector <- c(
    Edge_IND      = as.integer(E_true  >= lo_e_ind  & E_true  <= up_e_ind),
    Wedge_IND     = as.integer(S_true  >= lo_s_ind  & S_true  <= up_s_ind),
    Triangle_IND  = as.integer(T_true  >= lo_t_ind  & T_true  <= up_t_ind),
    CC_IND        = as.integer(CC_true >= lo_c_ind  & CC_true <= up_c_ind),
    Edge_EGO      = as.integer(E_true  >= lo_e_ego  & E_true  <= up_e_ego),
    Wedge_EGO     = as.integer(S_true  >= lo_s_ego  & S_true  <= up_s_ego),
    Triangle_EGO  = as.integer(T_true  >= lo_t_ego  & T_true  <= up_t_ego),
    CC_EGO        = as.integer(CC_true >= lo_c_ego  & CC_true <= up_c_ego)
  )
  
  results_table <- data.frame(
    Statistic = c("Edges","Wedges","Triangles","ClusteringCoeff"),
    True = c(E_true,S_true,T_true,CC_true),
    IND_Estimate = c(E_ind_est,S_ind_est,T_ind_est,C_ind_est),
    IND_Lower = c(lo_e_ind,lo_s_ind,lo_t_ind,lo_c_ind),
    IND_Upper = c(up_e_ind,up_s_ind,up_t_ind,up_c_ind),
    EGO_Estimate = c(E_ego_est,S_ego_est,T_ego_est,C_ego_est),
    EGO_Lower = c(lo_e_ego,lo_s_ego,lo_t_ego,lo_c_ego),
    EGO_Upper = c(up_e_ego,up_s_ego,up_t_ego,up_c_ego)
  )
  
  list(table = results_table,
       coverage = coverage_vector)
}

############################################################
# Run Evaluation
############################################################

results <- evaluate_sampling(A, A.ind, A.ego,
                             N, com.ind, Pi.mat.ind, p)

results
