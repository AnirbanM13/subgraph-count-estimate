rm(list = ls())

library(Rcpp)
library(RcppEigen)
library(knitr)

############################################################
## Simulation parameters
############################################################

M <- 2000
N.main <- 2 * 1e5

beta  <- 1/2
alpha <- 1/4

p.vec <- 0.04 / (N.main^alpha)

fname_sim_3 <- sprintf(
  "/sim-3-final/test-sim-3-iter-%d.Rdata",
  1:M
)

al.vec <- c(rep(1, N.main/2), rep(2, N.main/2))

Pi.mat <- rbind(
  c(1, 0.1),
  c(0.1, 3)
) / (N.main^beta)

############################################################
## Source variance file
############################################################

source("var-exp-final.R")

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
## Storage
############################################################

RES <- list(
  Edge = list(),
  Wedge = list(),
  Triangle = list(),
  CC = list()
)

iter_time <- numeric(M)
failed_iter <- list()

total_start <- Sys.time()

############################################################
## Main simulation loop 
############################################################

for (i in seq_along(fname_sim_3)) {
  
  iter_start <- Sys.time()
  
  tryCatch({
    
    dat <- get(load(fname_sim_3[i]))
    
    ########################################################
    ## EDGE
    ########################################################
    
    pop <- 2 * dat$all$edges
    est.ind <- 2 * dat$p1$edges1
    est.ego <- 2 * dat$p1$edges2
    
    ci.ind <- CI_eval(
      pop, est.ind,
      var_Edge_ind(N.main, al.vec, Pi.mat, p.vec),
      1 / p.vec^2
    )
    
    ci.ego <- CI_eval(
      pop, est.ego,
      var_Edge_ego(N.main, al.vec, Pi.mat, p.vec),
      1 / (2 * p.vec - p.vec^2)
    )
    
    RES$Edge[[length(RES$Edge) + 1]] <- data.frame(
      Method = c("Induced", "Ego-centric"),
      cover  = c(ci.ind$cover, ci.ego$cover),
      len    = c(ci.ind$len, ci.ego$len)
    )
    
    ########################################################
    ## WEDGE
    ########################################################
    
    pop <- 2 * dat$all$wedges
    est.ind <- 2 * dat$p1$wedges1
    est.ego <- 2 * dat$p1$wedges2
    
    ci.ind <- CI_eval(
      pop, est.ind,
      var_wedge_ind(N.main, al.vec, Pi.mat, p.vec),
      1 / p.vec^3
    )
    
    ci.ego <- CI_eval(
      pop, est.ego,
      var_wedge_ego(N.main, al.vec, Pi.mat, p.vec),
      1 / (p.vec + (1 - p.vec) * p.vec^2)
    )
    
    RES$Wedge[[length(RES$Wedge) + 1]] <- data.frame(
      Method = c("Induced", "Ego-centric"),
      cover  = c(ci.ind$cover, ci.ego$cover),
      len    = c(ci.ind$len, ci.ego$len)
    )
    
    ########################################################
    ## TRIANGLE
    ########################################################
    
    pop <- 6 * dat$all$triangles
    est.ind <- 6 * dat$p1$triangles1
    est.ego <- 6 * dat$p1$triangles2
    
    ci.ind <- CI_eval(
      pop, est.ind,
      var_Triangle_ind(N.main, al.vec, Pi.mat, p.vec),
      1 / p.vec^3
    )
    
    denom <- 3 * (1 - p.vec) * p.vec^2 + p.vec^3
    
    ci.ego <- CI_eval(
      pop, est.ego,
      var_Triangle_ego(N.main, al.vec, Pi.mat, p.vec),
      1 / denom
    )
    
    RES$Triangle[[length(RES$Triangle) + 1]] <- data.frame(
      Method = c("Induced", "Ego-centric"),
      cover  = c(ci.ind$cover, ci.ego$cover),
      len    = c(ci.ind$len, ci.ego$len)
    )
    
    ########################################################
    ## CLUSTERING COEFFICIENT
    ########################################################
    
    pop_cc <- CC_cnt(
      2 * dat$all$wedges,
      6 * dat$all$triangles
    )
    
    est.ind_cc <- CC_cnt(
      2 * dat$p1$wedges1,
      6 * dat$p1$triangles1
    )
    
    est.ego_cc <- CC_cnt(
      2 * dat$p1$wedges2,
      6 * dat$p1$triangles2
    )
    
    v.ind <- var_cc_ind(N.main, al.vec, Pi.mat, p.vec, beta)
    v.ego <- var_cc_ego(N.main, al.vec, Pi.mat, p.vec, beta)
    
    f <- (p.vec + (1 - p.vec) * p.vec^2) /
      (3 * p.vec^2 * (1 - p.vec) + p.vec^3)
    
    up.ind <- est.ind_cc + 1.96 * sqrt(v.ind) / (N.main^beta)
    lo.ind <- est.ind_cc - 1.96 * sqrt(v.ind) / (N.main^beta)
    
    up.ego <- f * est.ego_cc + 1.96 * sqrt(v.ego) / (N.main^beta)
    lo.ego <- f * est.ego_cc - 1.96 * sqrt(v.ego) / (N.main^beta)
    
    RES$CC[[length(RES$CC) + 1]] <- data.frame(
      Method = c("Induced", "Ego-centric"),
      cover = c(
        as.integer(lo.ind <= pop_cc & up.ind >= pop_cc),
        as.integer(lo.ego <= pop_cc & up.ego >= pop_cc)
      ),
      len = c(
        abs(up.ind - lo.ind),
        abs(up.ego - lo.ego)
      )
    )
    
    ########################################################
    ## Runtime
    ########################################################
    
    iter_time[i] <- as.numeric(
      difftime(Sys.time(), iter_start, units = "secs")
    )
    
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
## Total runtime
############################################################

total_time <- difftime(Sys.time(), total_start, units = "secs")
cat("\nTOTAL runtime:", round(total_time, 2), "seconds\n")

############################################################
## Summary builder
############################################################

make_table <- function(res, denom, stat) {
  
  R <- do.call(rbind, res)
  
  do.call(rbind, lapply(c("Induced", "Ego-centric"), function(m) {
    data.frame(
      Statistic = stat,
      Method = m,
      Result = paste0(
        mean(R$cover[R$Method == m]),
        " (",
        mean(R$len[R$Method == m]) / denom,
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
## Display + Save
############################################################

kable(
  FINAL_TABLE,
  caption = "Unified coverage and interval length results (Sim-3 sparse regime)"
)

write.table(
  FINAL_TABLE,
  "sim3_sparse_true_parameter_results.txt",
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)
