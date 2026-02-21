rm(list = ls())

library(Rcpp)
library(RcppEigen)
library(knitr)

############################################################
## Simulation parameters
############################################################

M <- 2000
N.main <- 2 * 1e6

alpha <- 1
beta  <- 0

c.val <- 50
p.val <- c.val / (N.main^alpha)

fname_sim_3 <- sprintf(
  "/home/anirbanm/Documents/sim-5-edge-0-1-final-1/test-%d.Rdata",
  1:M
)

al.vec <- c(rep(1, N.main/2), rep(2, N.main/2))

Pi.mat <- rbind(
  c(0.2, 0.05),
  c(0.05, 0.1)
) / (N.main^beta)

source("~/Desktop/r-file/var-exp-final.R")

############################################################
## Storage
############################################################

RES <- vector("list", M)
iter_time <- numeric(M)
failed_iter <- list()

total_start <- Sys.time()

############################################################
## Main simulation loop 
############################################################

for (i in seq_along(fname_sim_3)) {
  
  iter_start <- Sys.time()
  
  tryCatch({
    
    tmp <- get(load(fname_sim_3[i]))
    
    pop <- 2 * tmp$all$edges
    est.ego <- 2 * tmp$p1$edges2
    
    ########################################################
    ## Variance
    ########################################################
    
    v.ego <- var_Edge_ego(N.main, al.vec, Pi.mat, p.val)
    
    ########################################################
    ## Confidence Interval
    ########################################################
    
    denom <- 2 * p.val - p.val^2
    
    up <- (1 / denom) * (est.ego + 1.96 * sqrt(v.ego))
    lo <- (1 / denom) * (est.ego - 1.96 * sqrt(v.ego))
    
    ########################################################
    ## Coverage & Length
    ########################################################
    
    cover <- as.integer(lo <= pop & up >= pop)
    len   <- abs(up - lo)
    
    ########################################################
    ## Store
    ########################################################
    
    RES[[i]] <- data.frame(
      cover = cover,
      len   = len
    )
    
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
cat("Failed iterations:", length(failed_iter), "\n")

############################################################
## Final summary table
############################################################

ALL <- do.call(rbind, RES)

mean_cover <- mean(ALL$cover)
mean_len   <- mean(ALL$len) / (N.main^2)

FINAL_TABLE <- data.frame(
  Method = "Ego-centric",
  Result = paste0(
    round(mean_cover, 4),
    " (",
    mean_len,
    ")"
  )
)

############################################################
## Display + Save
############################################################

kable(
  FINAL_TABLE,
  caption = "Coverage and Interval Length (Edge, Ego only, p = c/N)"
)

write.table(
  FINAL_TABLE,
  "sim5_edge_ego_results.txt",
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)