rm(list = ls())

library(Rcpp)
library(RcppEigen)
library(knitr)

############################################################
## Simulation parameters
############################################################

M <- 2000
N.main <- 1e5
beta <- 1/2

fname_sim_2 <- sprintf(
  "/home/anirbanm/Documents/data-sim-2-final/test-sim-1-iter-%d.Rdata",
  1:M
)

############################################################
## Source external files
############################################################

source("~/Desktop/r-file/var-exp-final.R")
source("~/Desktop/r-file/bethe_hessian_clustering.R")
Rcpp::sourceCpp("Desktop/r-file/edge-list-to-adj-mat-conv-1.cpp")
Rcpp::sourceCpp("Desktop/r-file/pi_hat.cpp")

############################################################
## Helper functions
############################################################

CI_eval <- function(pop, est, var, scale = 1) {
  up <- scale * (est + 1.96 * sqrt(var))
  lo <- scale * (est - 1.96 * sqrt(var))
  list(
    cover = as.integer(lo <= pop & up >= pop),
    len   = abs(up - lo)
  )
}

CC_cnt <- function(w, t) ifelse(w > 0, t / w, 0)

############################################################
## Master storage
############################################################

RES <- list(
  Edge = list(),
  Wedge = list(),
  Triangle = list(),
  CC = list()
)

iter_time <- numeric(M)
failed_iter <- list()

############################################################
## TOTAL TIMER
############################################################

total_start <- Sys.time()

############################################################
## Main simulation loop (FAIL-SAFE)
############################################################

for (i in seq_along(fname_sim_2)) {
  
  iter_start <- Sys.time()
  
  tryCatch({
    
    dat <- get(load(fname_sim_2[i]))
    
    for (p in c(0.01, 0.05)) {
      
      tag <- ifelse(p == 0.01, "p1", "p2")
      
      adj <- dat[[tag]]$edge_list1[, 1:2]
      A   <- edge_list_to_adj_unique(adj)
      cl  <- bethe_hessian_clustering_adj(A)
      pi  <- PI_hat_cpp(cl, A)
      
      ########################################################
      ## EDGE
      ########################################################
      
      pop <- 2 * dat$all$edges
      est.ind <- 2 * dat[[tag]]$edges1
      est.ego <- 2 * dat[[tag]]$edges2
      
      ci.ind <- CI_eval(pop, est.ind,
                        var_Edge_ind(N.main, cl, pi, p),
                        1 / p^2)
      
      ci.ego <- CI_eval(pop, est.ego,
                        var_Edge_ego(N.main, cl, pi, p),
                        1 / (2 * p - p^2))
      
      RES$Edge[[length(RES$Edge) + 1]] <- data.frame(
        Method = c("Induced", "Ego-centric"),
        p = p,
        cover = c(ci.ind$cover, ci.ego$cover),
        len   = c(ci.ind$len, ci.ego$len)
      )
      
      ########################################################
      ## WEDGE
      ########################################################
      
      pop <- 2 * dat$all$wedges
      est.ind <- 2 * dat[[tag]]$wedges1
      est.ego <- 2 * dat[[tag]]$wedges2
      
      ci.ind <- CI_eval(pop, est.ind,
                        var_wedge_ind(N.main, cl, pi, p),
                        1 / p^3)
      
      ci.ego <- CI_eval(pop, est.ego,
                        var_wedge_ego(N.main, cl, pi, p),
                        1 / (p + (1 - p) * p^2))
      
      RES$Wedge[[length(RES$Wedge) + 1]] <- data.frame(
        Method = c("Induced", "Ego-centric"),
        p = p,
        cover = c(ci.ind$cover, ci.ego$cover),
        len   = c(ci.ind$len, ci.ego$len)
      )
      
      ########################################################
      ## TRIANGLE
      ########################################################
      
      pop <- 6 * dat$all$triangles
      est.ind <- 6 * dat[[tag]]$triangles1
      est.ego <- 6 * dat[[tag]]$triangles2
      
      ci.ind <- CI_eval(pop, est.ind,
                        var_Triangle_ind(N.main, cl, pi, p),
                        1 / p^3)
      
      ci.ego <- CI_eval(pop, est.ego,
                        var_Triangle_ego(N.main, cl, pi, p),
                        1 / (3 * (1 - p) * p^2 + p^3))
      
      RES$Triangle[[length(RES$Triangle) + 1]] <- data.frame(
        Method = c("Induced", "Ego-centric"),
        p = p,
        cover = c(ci.ind$cover, ci.ego$cover),
        len   = c(ci.ind$len, ci.ego$len)
      )
      
      ########################################################
      ## CLUSTERING COEFFICIENT (FIXED)
      ########################################################
      
      pop_cc <- CC_cnt(
        2 * dat$all$wedges,
        6 * dat$all$triangles
      )
      
      est.ind_cc <- CC_cnt(
        2 * dat[[tag]]$wedges1,
        6 * dat[[tag]]$triangles1
      )
      
      est.ego_cc <- CC_cnt(
        2 * dat[[tag]]$wedges2,
        6 * dat[[tag]]$triangles2
      )
      
      v.ind <- var_cc_ind(N.main, cl, pi, p, beta)
      v.ego <- var_cc_ego(N.main, cl, pi, p, beta)
      
      f <- (p + (1 - p) * p^2) / (3 * p^2 * (1 - p) + p^3)
      
      ## Induced CI
      up.ind <- est.ind_cc + 1.96 * sqrt(v.ind) / (N.main^beta)
      lo.ind <- est.ind_cc - 1.96 * sqrt(v.ind) / (N.main^beta)
      
      ## Ego CI (bias-corrected)
      up.ego <- f * est.ego_cc + 1.96 * sqrt(v.ego) / (N.main^beta)
      lo.ego <- f * est.ego_cc - 1.96 * sqrt(v.ego) / (N.main^beta)
      
      RES$CC[[length(RES$CC) + 1]] <- data.frame(
        Method = c("Induced", "Ego-centric"),
        p = p,
        cover = c(
          as.integer(lo.ind <= pop_cc & up.ind >= pop_cc),
          as.integer(lo.ego <= pop_cc & up.ego >= pop_cc)
        ),
        len = c(
          abs(up.ind - lo.ind),
          abs(up.ego - lo.ego)
        )
      )
    }
    
    iter_time[i] <- as.numeric(difftime(Sys.time(), iter_start, units = "secs"))
    cat("Iteration", i, ":", round(iter_time[i], 2), "seconds\n")
    
  }, error = function(e) {
    
    iter_time[i] <- NA
    failed_iter[[length(failed_iter) + 1]] <<- list(
      iteration = i,
      error = conditionMessage(e)
    )
    
    cat("Iteration", i, ": FAILED\n")
  })
}

############################################################
## TOTAL TIME
############################################################

total_time <- difftime(Sys.time(), total_start, units = "secs")
cat("\nTOTAL runtime:", round(total_time, 2), "seconds\n")
cat("Failed iterations:", length(failed_iter), "\n")


############################################################
## Save runtime + failure logs
############################################################

write.table(
  data.frame(
    iteration = 1:M,
    runtime_seconds = iter_time
  ),
  "iteration_runtime_log.txt",
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

if (length(failed_iter) > 0) {
  write.table(
    do.call(rbind, lapply(failed_iter, as.data.frame)),
    "failed_iterations_log.txt",
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
}

############################################################
## Summary table builder
############################################################

make_table <- function(res, denom = 1, stat) {
  
  R <- do.call(rbind, res)
  
  do.call(rbind, lapply(c("Induced", "Ego-centric"), function(m) {
    data.frame(
      Statistic = stat,
      Method = m,
      `p = 0.01` = paste0(
        mean(R$cover[R$Method == m & R$p == 0.01]),
        " (",
        mean(R$len[R$Method == m & R$p == 0.01]) / denom,
        ")"
      ),
      `p = 0.05` = paste0(
        mean(R$cover[R$Method == m & R$p == 0.05]),
        " (",
        mean(R$len[R$Method == m & R$p == 0.05]) / denom,
        ")"
      )
    )
  }))
}

############################################################
## Final unified table
############################################################

FINAL_TABLE <- rbind(
  make_table(RES$Edge,     N.main^2, "Edge"),
  make_table(RES$Wedge,    N.main^3, "Wedge"),
  make_table(RES$Triangle, N.main^3, "Triangle"),
  make_table(RES$CC,       1,        "Clustering Coefficient")
)

############################################################
## Display + save
############################################################

kable(
  FINAL_TABLE,
  caption = "Unified coverage and interval length results (Sim-2, all statistics,Parameter unknown)"
)

write.table(
  FINAL_TABLE,
  "sim2_all_statistics_est_1_results.txt",
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)
