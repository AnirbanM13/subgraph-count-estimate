rm(list= ls())



library(Rcpp)


sourceCpp("sim-5-data-gen-engine.cpp")


N.main <- 2*1e6
al.vec <- c(rep(1,N.main/2), rep(2,N.main/2))
alpha = 1
beta = 0
Pi.mat <- rbind(c(0.2, 0.05), c(0.05, 0.1))/N.main^(beta)
p.vec = 50/(N.main^(alpha))
#p.vec <- c(0.01,0.05)
M = 2000

# Print start time
overall_start <- Sys.time()
cat("Code started at:", format(overall_start), "\n")

for (i in 1:M) {
  t0 <- Sys.time()
  
  list.full <- generate_sbm_edges_fast_no_edgelist(
    al.vec, Pi.mat, p.vec,
    threads = 64,
    seed = as.integer(Sys.time()) + i
  )
  
  fname.i <- paste("test-", i, ".Rdata", sep = "")
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
