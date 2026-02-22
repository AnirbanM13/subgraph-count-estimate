rm(list=ls())
library(Rcpp)

Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
Sys.setenv("PKG_LIBS"="-fopenmp")
sourceCpp("sim-1-2-3-data-gen-engine.cpp")

N.main <- 1e5
al.vec <- c(rep(1,N.main/2), rep(2,N.main/2))
beta = 1/2
Pi.mat <- rbind(c(1.6, 0.65), c(0.65, 2.85))/(N.main^(beta))
#Pi.mat <- rbind(c(0.2, 0.05), c(0.05, 0.1))

p.vec = c(1, 5)/100
M = 2000

# Print start time
overall_start <- Sys.time()
cat("Code started at:", format(overall_start), "\n")

for (i in 1:M) {
  t0 <- Sys.time()
  
  list.full <- generate_sbm_edges_with_multiple_W_fast(
    al.vec, Pi.mat, p.vec,
    threads = 64,
    seed = as.integer(Sys.time()) + i
  )
  
  fname.i <- paste("test-sim-1-iter-", i, ".Rdata", sep = "")
  save(list.full, file = fname.i)
  
  t1 <- Sys.time()
  elapsed <- difftime(t1, t0, units = "secs")
  
  cat("Iteration m =", i, "done in", elapsed, "seconds (finished at", format(t1), ")\n")
}

# Print end time
overall_end <- Sys.time()
total_elapsed <- difftime(overall_end, overall_start, units = "secs")
cat("Code ended at:", format(overall_end), "\n")
cat("Total runtime:", total_elapsed, "seconds\n")


