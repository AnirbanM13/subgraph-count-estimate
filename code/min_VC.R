# Function to check if a set of vertices is a vertex cover
is_vertex_cover <- function(adj_matrix, vertices) {
  n <- nrow(adj_matrix)
  
  for (i in 1:n) {
    for (j in 1:n) {
      if (adj_matrix[i, j] == 1) {
        # if edge exists but neither endpoint is in the set
        if (!(i %in% vertices || j %in% vertices)) {
          return(FALSE)
        }
      }
    }
  }
  return(TRUE)
}

# Function to find all minimum vertex covers
find_min_vertex_covers <- function(adj_matrix) {
  n <- nrow(adj_matrix)
  vertices <- 1:n
  
  min_size <- n
  min_covers <- list()
  
  # generate all subsets
  for (k in 1:n) {
    subsets <- combn(vertices, k, simplify = FALSE)
    
    for (subset in subsets) {
      if (is_vertex_cover(adj_matrix, subset)) {
        
        if (k < min_size) {
          min_size <- k
          min_covers <- list(subset)
        } else if (k == min_size) {
          min_covers <- append(min_covers, list(subset))
        }
      }
    }
    
    # stop early once minimum size is found
    if (length(min_covers) > 0) {
      break
    }
  }
  
  return(min_covers)
}

# Example adjacency matrix
adj <- matrix(c(
  0,1,1,0,0,
  1,0,1,0,0,
  1,1,0,1,1,
  0,0,1,0,1,
  0,0,1,1,0
), nrow = 5, byrow = TRUE)

# Run
covers <- find_min_vertex_covers(adj1)

# Print results
print("Minimum Vertex Covers:")
print(covers)
