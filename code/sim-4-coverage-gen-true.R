rm(list = ls())

library(Rcpp)
library(RcppEigen)
library(knitr)

############################################################
## Simulation parameters
############################################################

M <- 2000
N.main <- 1e6

alpha <- 0
beta  <- 2

p.vec <- c(0.01, 0.05)

fname_sim_3 <- sprintf(
  "/home/anirbanm/Documents/sim-edge-0-2-final/test-%d.Rdata",
  1:M
)

al.vec <- c(rep(1, N.main/2), rep(2, N.main/2))

Pi.mat <- rbind(
  c(3e5, 1e5),
  c(1e5, 2e5)
) / (N.main^beta)

############################################################
## Master storage (ONLY FOR EDGE)
############################################################

RES <- list()
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
    pop <- tmp$all$edges
    
    ############################################
    ## Helper for one p
    ############################################
    
    run_one_p <- function(p, est.ind, est.ego) {
      
      ########################
      ## Induced
      ########################
      
      poi_1 <- 0.5 * 0.25 * p^2 * sum(N.main^beta * Pi.mat)
      poi_2 <- 0.5 * 0.25 * (1 - p^2) * sum(N.main^beta * Pi.mat)
      
      U1 <- rpois(1e6, poi_1)
      U2 <- rpois(1e6, poi_2)
      
      Z  <- (1 - p^2) * U1 - p^2 * U2
      
      lower.ind <- quantile(Z, 0.025, names = FALSE)
      upper.ind <- quantile(Z, 0.975, names = FALSE)
      
      up.ind <- (est.ind - lower.ind) / p^2
      lo.ind <- (est.ind - upper.ind) / p^2
      
      ########################
      ## Ego-centric
      ########################
      
      d <- 2 * p - p^2
      
      poi_1e <- 0.5 * 0.25 * d * sum(N.main^beta * Pi.mat)
      poi_2e <- 0.5 * 0.25 * (1 - d) * sum(N.main^beta * Pi.mat)
      
      U1e <- rpois(1e6, poi_1e)
      U2e <- rpois(1e6, poi_2e)
      
      Ze  <- (1 - d) * U1e - d * U2e
      
      lower.ego <- quantile(Ze, 0.025, names = FALSE)
      upper.ego <- quantile(Ze, 0.975, names = FALSE)
      
      up.ego <- (est.ego - lower.ego) / d
      lo.ego <- (est.ego - upper.ego) / d
      
      ########################
      ## Coverage & Length
      ########################
      
      list(
        cover.ind = as.integer(lo.ind <= pop & up.ind >= pop),
        cover.ego = as.integer(lo.ego <= pop & up.ego >= pop),
        len.ind   = abs(up.ind - lo.ind),
        len.ego   = abs(up.ego - lo.ego)
      )
    }
    
    ########################################################
    ## p = 0.01
    ########################################################
    
    out1 <- run_one_p(
      p.vec[1],
      tmp$p1$edges1,
      tmp$p1$edges2
    )
    
    ########################################################
    ## p = 0.05
    ########################################################
    
    out2 <- run_one_p(
      p.vec[2],
      tmp$p2$edges1,
      tmp$p2$edges2
    )
    
    ########################################################
    ## Store
    ########################################################
    
    RES[[i]] <- data.frame(
      p = rep(p.vec, each = 2),
      Method = rep(c("Induced","Ego-centric"), 2),
      cover = c(
        out1$cover.ind,
        out1$cover.ego,
        out2$cover.ind,
        out2$cover.ego
      ),
      len = c(
        out1$len.ind,
        out1$len.ego,
        out2$len.ind,
        out2$len.ego
      )
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

############################################################
## Summary table builder
############################################################

ALL <- do.call(rbind, RES)

FINAL_TABLE <- do.call(rbind, lapply(c("Induced","Ego-centric"), function(m) {
  
  data.frame(
    Method = m,
    `p = 0.01` = paste0(
      mean(ALL$cover[ALL$Method==m & ALL$p==0.01]),
      " (",
      mean(ALL$len[ALL$Method==m & ALL$p==0.01]) / N.main^2,
      ")"
    ),
    `p = 0.05` = paste0(
      mean(ALL$cover[ALL$Method==m & ALL$p==0.05]),
      " (",
      mean(ALL$len[ALL$Method==m & ALL$p==0.05]) / N.main^2,
      ")"
    )
  )
}))

############################################################
## Display + Save
############################################################

kable(
  FINAL_TABLE,
  caption = "Coverage and Interval Length (Edge, Poisson CI, N=1e6)"
)

write.table(
  FINAL_TABLE,
  "sim4_edge_poisson_true_parameter_results.txt",
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)