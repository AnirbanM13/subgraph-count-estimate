construct_G_t_H <- function(adj_matrix, t) {
  n <- nrow(adj_matrix)
  
  # Step 1: generate all t-subsets (tuples jbb)
  combos <- combn(n, t, simplify = FALSE)
  
  subgraphs <- list()
  edge_counts <- c()
  
  # Step 2: compute induced subgraphs and edge counts
  for (i in seq_along(combos)) {
    v <- combos[[i]]
    
    sub_adj <- adj_matrix[v, v]
    edges <- sum(sub_adj) / 2   # undirected graph
    
    subgraphs[[i]] <- list(
      vertices = v,
      adj = sub_adj,
      edges = edges
    )
    
    edge_counts[i] <- edges
  }
  
  # Step 3: compute m(t:H)
  max_edges <- max(edge_counts)
  
  # Step 4: select extremal ones
  G_t_H <- subgraphs[edge_counts == max_edges]
  
  return(list(
    max_edges = max_edges,
    G_t_H = G_t_H
  ))
}


adj <- matrix(c(
  0,1,1,0,0,
  1,0,1,0,0,
  1,1,0,1,1,
  0,0,1,0,0,
  0,0,1,0,0
), nrow = 5, byrow = TRUE)



adj1 <- matrix(c(
  0,1,1,0,0,0,
  1,0,1,0,0,0,
  1,1,0,1,0,0,
  0,0,1,0,1,1,
  0,0,0,1,0,1,
  0,0,0,1,1,0
), nrow = 6, byrow = TRUE)

result <- construct_G_t_H(adj1, t = 4)

print(result$max_edges)
print(result$G_t_H)
