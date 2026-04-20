# Function needed for checking assumption 5
## Overview

This project provides a framework for analyzing structural properties of graphs using adjacency matrices. It focuses on:

* Constructing extremal induced subgraphs ( G_t(H) )
* Computing minimum vertex covers
* Studying inclusion relationships between vertex sets
* Verifying the **A5 condition** across all subset sizes

---

## Core Idea

Given a graph ( H ) and an integer ( t ):

1. Generate all subsets of vertices of size ( t )
2. Find induced subgraphs on these subsets
3. Select those with the **maximum number of edges**
4. Compare these vertex sets with **minimum vertex covers**

---

## Functions

---

### 1. `construct_G_t_H(adj_matrix, t)`

Constructs the extremal induced subgraphs ( G_t(H) ).

#### Steps

1. Generate all ( t )-vertex subsets
2. Compute induced subgraphs
3. Count edges in each subgraph
4. Select subgraphs with maximum edge count

#### Returns

* `max_edges`: Maximum number of edges among all ( t )-subgraphs
* `G_t_H`: List of extremal subgraphs, each containing:

  * `vertices`: Vertex indices
  * `adj`: Adjacency matrix of subgraph
  * `edges`: Number of edges

#### Example

```r
result <- construct_G_t_H(adj1, t = 4)

print(result$max_edges)
print(result$G_t_H)
```

---

### 2. `analyze_graph_relation(adj_matrix, t)`

Analyzes relationships between:

* Extremal vertex sets ( V \in G_t(H) )
* Minimum vertex covers ( D )

#### Cases

* **If ( t > \tau(H) )** → Check: ( D \subseteq V )
* **If ( t < \tau(H) )** → Check: ( V \subseteq D )
* **If ( t = \tau(H) )** → Check: ( D = V )

#### Returns

* Case type
* Minimum vertex cover size ( \tau(H) )
* All extremal vertex sets
* All minimum vertex covers
* Matching pairs satisfying conditions

---

### 3. `check_A5_condition(adj_matrix)`

Checks whether graph ( H ) satisfies the **A5 condition**.

#### Logic

* Iterate over all ( t = 1 ) to ( n )
* Run `analyze_graph_relation`
* If any ( t ) fails → A5 is **not satisfied**

#### Output

```r
A5 is satisfied for all t
```

or

```r
A5 is NOT satisfied (failed at t = X)
```

---

### 4. `path_adj_matrix(n)`

Generates adjacency matrix of path graph ( P_n ).

---

## Example Graphs

### Graph 1 (5 vertices)

```r
adj <- matrix(c(
  0,1,1,0,0,
  1,0,1,0,0,
  1,1,0,1,1,
  0,0,1,0,0,
  0,0,1,0,0
), nrow = 5, byrow = TRUE)
```

### Graph 2 (6 vertices)

```r
adj1 <- matrix(c(
  0,1,1,0,0,0,
  1,0,1,0,0,0,
  1,1,0,1,0,0,
  0,0,1,0,1,1,
  0,0,0,1,0,1,
  0,0,0,1,1,0
), nrow = 6, byrow = TRUE)
```

---

## Example Usage

```r
# Construct G_t(H)
result <- construct_G_t_H(adj1, t = 4)

# Analyze relations
analyze_graph_relation(adj1, t = 4)

# Check A5 condition
check_A5_condition(adj1)
```

---

## Dependencies

You must provide:

* `find_min_vertex_covers(adj_matrix)`

---

## Complexity Notes

* Generating all subsets: ( O(\binom{n}{t}) )
* Vertex cover computation: NP-hard
* Overall runtime grows rapidly with graph size

---

## Interpretation of Results

| Relation | Meaning                            |
| -------- | ---------------------------------- |
| D ⊆ V    | Cover is contained in extremal set |
| V ⊆ D    | Extremal set lies inside cover     |
| D = V    | Exact match                        |

---

## Summary

This project connects:

* Extremal subgraph theory
* Vertex cover structure
* Combinatorial conditions (A5)

It can be used for:

* Theoretical graph research
* Testing conjectures
* Structural graph analysis

---


# Assumption 5 Checker

## Overview

This project provides tools for analyzing structural relationships in graphs using their adjacency matrices. It focuses on comparing:

* Vertex sets derived from a constructed graph ( G_t(H) )
* Minimum vertex covers of the original graph ( H )

It also verifies whether a graph satisfies a specific theoretical condition referred to as **A5**.

---

## Features

* Construct and analyze ( G_t(H) ) (external dependency)
* Compute minimum vertex covers of a graph
* Compare relationships between:

  * Vertex sets ( V )
  * Minimum vertex covers ( D )
* Automatically classify cases based on:

  * ( t > \tau(H) )
  * ( t < \tau(H) )
  * ( t = \tau(H) )
* Validate the **A5 condition** across all values of ( t )

---

## Definitions

* **Adjacency Matrix**: Matrix representation of a graph.
* **Vertex Cover**: A set of vertices such that every edge is incident to at least one vertex in the set.
* **Minimum Vertex Cover**: Vertex cover of smallest possible size.
* **( \tau(H) )**: Size of a minimum vertex cover of graph ( H ).
* **( G_t(H) )**: A derived graph depending on parameter ( t ) (constructed via `construct_G_t_H()`).

---

## Functions

### 1. `analyze_graph_relation(adj_matrix, t)`

Analyzes the relationship between:

* Vertex sets from ( G_t(H) )
* Minimum vertex covers of ( H )

#### Behavior

Depending on ( t ) and ( \tau(H) ):

* **Case 1:** ( t > \tau(H) ) → checks if ( D \subseteq V )
* **Case 2:** ( t < \tau(H) ) → checks if ( V \subseteq D )
* **Case 3:** ( t = \tau(H) ) → checks if ( D = V )

#### Returns

A list containing:

* `type`: Case classification
* `tau`: Minimum vertex cover size
* `Gt_vertices`: Vertex sets from ( G_t(H) )
* `min_vertex_covers`: All minimum vertex covers
* `matches`: Valid pairs satisfying the condition

---

### 2. `check_A5_condition(adj_matrix)`

Checks whether the graph satisfies the **A5 condition**.

#### Logic

* Iterates over all ( t = 1 ) to ( n ) (number of vertices)
* Calls `analyze_graph_relation()` for each ( t )
* If any ( t ) fails → A5 is **not satisfied**
* If all pass → A5 is **satisfied**

#### Returns

* `TRUE` → A5 satisfied
* `FALSE` → A5 not satisfied

---

### 3. `path_adj_matrix(n)`

Generates the adjacency matrix for a **path graph ( P_n )**.

#### Example

```
P4 <- path_adj_matrix(4)
```

---

## Example Usage

```r
# Create path graphs
P4 <- path_adj_matrix(4)
P5 <- path_adj_matrix(5)
P6 <- path_adj_matrix(6)

# Check A5 condition
check_A5_condition(P4)
check_A5_condition(P5)
check_A5_condition(P6)
```

---

## Dependencies

The following functions must be implemented separately:

* `construct_G_t_H(adj_matrix, t)`
* `find_min_vertex_covers(adj_matrix)`

---

## Output Interpretation

* `"D ⊆ V"` → Minimum vertex cover is contained in ( V )
* `"V ⊆ D"` → Vertex set is contained in minimum cover
* `"D = V"` → Exact match

If no matches are found:

```
No V satisfies the condition.
```

---

## Notes

* Performance depends heavily on the implementation of:

  * Vertex cover computation
  * ( G_t(H) ) construction
* Minimum vertex cover is an NP-hard problem, so large graphs may be computationally expensive.



# var_H: Variance Computation 

## Overview

This project implements an efficient computation of the variance term `var_H` for graph-based statistical models. It combines:

* **C++ (via Rcpp)** for performance-critical combinatorial and summation routines
* **R** for orchestration, graph processing, and higher-level statistical logic
* **Parallel computing** to handle combinatorial explosion

The implementation is designed for **high-dimensional graph structures** and supports both:

* `"ego-centric"` mode (default)
* `"induced"` mode

---

## Key Components

### 1. `generate_Rt_cpp`

Generates all possible tuples `(j, k, xi)` where:

* `j`, `k`: combinations of size `t` from `{1,...,R}`
* `xi`: permutations of `{1,...,t}`

This defines the **indexing structure** for summation terms.

---

### 2. `theta_model_sum_cpp` (B-term)

A C++-accelerated function that computes:

* A **large combinatorial sum** over assignments of latent variables `u` and `v`
* Uses **depth-first search (DFS)** to enumerate all configurations
* Applies:

  * Edge-based probability weights (`Pi`)
  * Node weights (`lambda`)
  * Structural constraints from `(j, k, xi)`

This is the **core probabilistic computation**.

---

### 3. Vertex Cover Functions (A-term)

These functions operate on graph structures:

#### `H1_union_H2`

* Builds a merged graph using union-find (disjoint set)
* Constructs a compressed graph representation
* Computes **all vertex covers** of the resulting graph

#### `psi_design_cov`

* Computes a weighted sum over vertex covers:

  $
  \sum p^{|C|}(1-p)^{(2R - t - |C|)}
  $

#### `get_all_vertex_covers`

* Enumerates all vertex covers of a graph via bitmasking

#### `f2_vertex_cover`

* Computes a baseline vertex cover probability term

---

### 4. `var_H` (Main Function)

This is the final function that computes the variance.

### Function Signature

```r
var_H(
  adj_matrix,  # Adjacency matrix of graph H
  Pi,          # Edge probability matrix (K x K)
  lambda,      # Node distribution vector (length K)
  N,           # Network size scaling parameter
  pN,          # Sampling/design probability
  mode = c("ego-centric", "induced"),
  mc.cores = 8 # Number of cores for parallelization
)
```

---

## How It Works

For each `t = 1,...,R`:

1. Generate all `(j, k, xi)` using `generate_Rt_cpp`
2. For each triple:

   * Compute **B-term** using `theta_model_sum_cpp`
   * If `"ego-centric"`:

     * Compute graph union and vertex covers
     * Compute **A-term correction**
     * Multiply `(A - baseline²) * B`
3. Aggregate results across all `t`
4. Apply scaling based on:

   * `N` (network size)
   * `pN` (sampling probability)
   * Mode selection

---

## Modes

### `"ego-centric"` (default)

* Includes both **A-term (design correction)** and **B-term**
* Captures sampling variability and graph structure

### `"induced"`

* Only computes the **B-term**
* Faster but ignores design-based corrections

---

## Performance Considerations

* The algorithm is **computationally expensive**:

  * Combination count grows as ( O(\binom{R}{t}^2 \cdot t!) )
  * DFS over assignments is exponential in `R`
* Parallelization via `mclapply` is essential
* C++ (Rcpp) is used to handle inner loops efficiently

---

## Example Usage

```r
tri <- matrix(
  c(0,1,1,
    1,0,1,
    1,1,0),
  nrow = 3, byrow = TRUE
)

Pi.mat <- rbind(c(0.2, 0.05),
                c(0.05, 0.1))

lambda <- c(0.5, 0.5)

result <- var_H(
  adj_matrix = tri,
  Pi = Pi.mat,
  lambda = lambda,
  N = 20000,
  pN = 0.01,
  mode = "ego-centric",
  mc.cores = 4
)
```

---

## Dependencies

* R packages:

  * `Rcpp`
  * `parallel`

---

## Notes

* Suitable for **small to moderate graph sizes** (due to exponential complexity)
* For large graphs, approximations or sampling-based methods may be needed
* The implementation emphasizes **exact computation over approximation**

---

## Summary

`var_H` computes a **variance functional over graph structures** by combining:

* Combinatorics (combinations & permutations)
* Graph theory (vertex covers, unions)
* Probabilistic modeling (latent assignments, edge probabilities)
* High-performance computing (Rcpp + parallelization)

---

## Author Notes

This implementation is designed for **research-level statistical modeling** where exact combinatorial evaluation is required.

---
