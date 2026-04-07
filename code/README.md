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

  [
  \sum p^{|C|}(1-p)^{(2R - t - |C|)}
  ]

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
