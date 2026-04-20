analyze_graph_relation <- function(adj_matrix, t) {
  
  # ---- Step 1: Compute G_t(H) ----
  Gt <- construct_G_t_H(adj_matrix, t)
  Gt_vertices <- lapply(Gt$G_t_H, function(x) x$vertices)
  
  # ---- Step 2: Minimum vertex covers ----
  min_covers <- find_min_vertex_covers(adj_matrix)
  tau <- length(min_covers[[1]])
  
  results <- list()
  
  # ---- Case 1: t > tau(H) ----
  if (t > tau) {
    case_type <- "t > tau(H): D ⊆ V"
    
    for (V in Gt_vertices) {
      for (D in min_covers) {
        if (all(D %in% V)) {
          results <- append(results, list(list(
            V = V,
            D = D,
            relation = "D ⊆ V"
          )))
          break
        }
      }
    }
  }
  
  # ---- Case 2: t < tau(H) ----
  else if (t < tau) {
    case_type <- "t < tau(H): V ⊆ D"
    
    for (V in Gt_vertices) {
      for (D in min_covers) {
        if (all(V %in% D)) {
          results <- append(results, list(list(
            V = V,
            D = D,
            relation = "V ⊆ D"
          )))
          break
        }
      }
    }
  }
  
  # ---- Case 3: t = tau(H) ----
  else {
    case_type <- "t = tau(H): D = V"
    
    for (V in Gt_vertices) {
      for (D in min_covers) {
        if (setequal(V, D)) {
          results <- append(results, list(list(
            V = V,
            D = D,
            relation = "D = V"
          )))
        }
      }
    }
  }
  
  if (length(results) == 0) {
    cat("No V satisfies the condition.\n")
  }
  
  return(list(
    type = case_type,
    tau = tau,
    Gt_vertices = Gt_vertices,
    min_vertex_covers = min_covers,
    matches = results
  ))
}




# result <- analyze_graph_relation(adj1, t = 5)
# 
# print(result$type)
# print(result$matches)
# print(result$tau)




check_A5_condition <- function(adj_matrix) {
  
  n <- nrow(adj_matrix)
  
  for (t in 1:n) {
    
    cat("Checking t =", t, "...\n")
    
    result <- analyze_graph_relation(adj_matrix, t)
    
    if (length(result$matches) == 0) {
      cat("A5 is NOT satisfied (failed at t =", t, ")\n")
      return(FALSE)
    }
  }
  
  cat("A5 is satisfied for all t\n")
  return(TRUE)
}






#Example 
# Function to create adjacency matrix of path graph P_n
path_adj_matrix <- function(n) {
  A <- matrix(0, n, n)
  for (i in 1:(n-1)) {
    A[i, i+1] <- 1
    A[i+1, i] <- 1
  }
  return(A)
}

# Matrices
P4 <- path_adj_matrix(4)
P5 <- path_adj_matrix(5)
P6 <- path_adj_matrix(6)

# Print
P4
P5
P6




check_A5_condition(P4)
check_A5_condition(P5)
check_A5_condition(P6)


