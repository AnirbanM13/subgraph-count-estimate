rm(list = ls())

library(Rcpp)
library(RcppEigen)
library(knitr)

############################################################
## Simulation parameters
############################################################

M <- 2000
N.main <- 20000

fname_sim_1 <- sprintf(
  "/data-sim-1-final/sim-1-iter-%d.Rdata",
  1:M
)

al.vec <- c(rep(1, N.main/2), rep(2, N.main/2))
Pi.mat <- rbind(c(0.2, 0.05), c(0.05, 0.1))

############################################################
## Source files
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

############################################################
## TOTAL TIMER
############################################################

total_start <- Sys.time()

############################################################
## Main simulation loop
############################################################

for (i in seq_along(fname_sim_1)) {
  
  iter_start <- Sys.time()
  
  tryCatch({
    
    dat <- get(load(fname_sim_1[i]))
    
    for (p in c(0.01, 0.05)) {
      
      tag <- ifelse(p == 0.01, "p1", "p2")
      
      ########################################################
      ## EDGE
      ########################################################
      
      pop <- 2 * dat$all$edges
      est.ind <- 2 * dat[[tag]]$edges1
      est.ego <- 2 * dat[[tag]]$edges2
      
      ci.ind <- CI_eval(pop, est.ind,
                        var_Edge_ind(N.main, al.vec, Pi.mat, p),
                        1 / p^2)
      
      ci.ego <- CI_eval(pop, est.ego,
                        var_Edge_ego(N.main, al.vec, Pi.mat, p),
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
                        var_wedge_ind(N.main, al.vec, Pi.mat, p),
                        1 / p^3)
      
      ci.ego <- CI_eval(pop, est.ego,
                        var_wedge_ego(N.main, al.vec, Pi.mat, p),
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
                        var_Triangle_ind(N.main, al.vec, Pi.mat, p),
                        1 / p^3)
      
      ci.ego <- CI_eval(pop, est.ego,
                        var_Triangle_ego(N.main, al.vec, Pi.mat, p),
                        1 / (3 * (1 - p) * p^2 + p^3))
      
      RES$Triangle[[length(RES$Triangle) + 1]] <- data.frame(
        Method = c("Induced", "Ego-centric"),
        p = p,
        cover = c(ci.ind$cover, ci.ego$cover),
        len   = c(ci.ind$len, ci.ego$len)
      )
      
      ########################################################
      ## CLUSTERING COEFFICIENT 
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
      
      v.ind <- var_cc_ind(N.main, al.vec, Pi.mat, p, beta = 0)
      v.ego <- var_cc_ego(N.main, al.vec, Pi.mat, p, beta = 0)
      
      f <- (p + (1 - p) * p^2) / (3 * p^2 * (1 - p) + p^3)
      
      up.ind <- est.ind_cc + 1.96 * sqrt(v.ind)
      lo.ind <- est.ind_cc - 1.96 * sqrt(v.ind)
      
      up.ego <- f * est.ego_cc + 1.96 * sqrt(v.ego)
      lo.ego <- f * est.ego_cc - 1.96 * sqrt(v.ego)
      
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
## Summary builder
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
## Display + Save
############################################################

kable(
  FINAL_TABLE,
  caption = "Unified coverage and interval length results (Sim-1, true parameters)"
)

write.table(
  FINAL_TABLE,
  "sim1_true_parameter_results.txt",
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)
