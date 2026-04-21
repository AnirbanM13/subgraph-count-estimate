# Prediction Intervals for Subgraph Densities in Sampled Networks

## Overview

This project studies the construction and empirical evaluation of **asymptotic prediction intervals (PIs)** for:

* Subgraph densities (edge, wedge, triangle)
* Clustering coefficient

within large random networks modeled using the **Stochastic Block Model (SBM)** under various **sampling schemes** and **sparsity regimes**.

---

## Prediction Interval Definition

Let:

* $T_{N,l}(H) \Rightarrow F$, where $F$ is a cumulative distribution function (c.d.f.)
* $H$: subgraph of interest
* $R = |V(H)|$: number of vertices in $H$

Then, an asymptotic two-sided level $(1 - \eta)$ prediction interval for the population subgraph density $N^{-R} \cdot S_N(H)$ is:

$$
J_{\eta,l,F}(H) \equiv
\left[
\frac{\widehat{S}*{N,l}(H) - \xi*{1-\eta/2} \cdot \sigma_{N,l}(H)}{N^R \cdot f_l(p_N:H)},
\frac{\widehat{S}*{N,l}(H) - \xi*{\eta/2} \cdot \sigma_{N,l}(H)}{N^R \cdot f_l(p_N:H)}
\right]
$$

Where:

* $\xi_\eta$: $\eta$-quantile of distribution $F$
* $f_l(p_N:H)$: sampling correction factor
* $\sigma_{N,l}(H)$: asymptotic standard deviation

---

## Objectives

* Evaluate **empirical coverage** of prediction intervals
* Compare performance under:

  * Different **subgraphs**: edge, wedge, triangle
  * Different **sampling methods**:

    * Induced sampling
    * Ego-centric sampling
  * Various **network sparsity levels**

---

## Key Challenge

The main difficulty lies in estimating:

$$
\sigma^2_{N,l}(H)
$$

This depends on **unknown model parameters**, including:

* Limiting matrix $\mathbf{C}$
* Class proportions $(\lambda_1, \ldots, \lambda_K)$
* Edge probabilities $\pi_{N,u,v}$
* Sparsity parameter $\beta$

### Assumptions

* Model sparsity $\beta$ is assumed **known**
* Sampling sparsity $\alpha$ is assumed **known**

Misestimating $\beta$ may lead to choosing the wrong limiting distribution (Gaussian vs Poisson), affecting PI accuracy.

---

## Parameter Estimation

### Community Detection

We use the **Bethe-Hessian clustering algorithm** for detecting communities:

* Estimates number of communities $\widehat{K}$
* Assigns class labels to sampled nodes

---

### Edge Probability Estimation

$$
\widehat{\pi}*{k,m}(l) =
\frac{
\sum*{(i,j)\in I_{l,N}} Y_{i,j} \cdot \mathbf{1}(\widehat{\alpha}*i=k,\widehat{\alpha}*j=m)
}{
\sum*{(i,j)\in I*{l,N}} \mathbf{1}(\widehat{\alpha}_i=k,\widehat{\alpha}_j=m)
}
$$

#### Notes:

* $l = 1$: Induced sampling
* $l = 2$: Ego-centric sampling

Ego-centric case uses:

* Labels from induced data
* Edge estimates from ego-centric data

---

## Theoretical Limitation

Due to **dependence in sampled data**, standard guarantees (e.g., consistency) may fail.

As a result:

> The estimated prediction interval $\widehat{J}_{\eta,l,F}(H)$ has no theoretical guarantee of correct coverage.

---

## Simulation Study

### Setup

* **2000 Monte Carlo replications**

Each replication:

1. Generate full network $\mathbf{Y}_N$
2. Sample nodes via Bernoulli sampling
3. Observe sampled network (induced / ego-centric)

---

### Common Configuration

* $K = 2$ communities
* Equal class proportions:

$$
\lambda_1 = \lambda_2 = \frac{1}{2}
$$

---

## Simulation Scenarios

### (S.1) Dense Network & Dense Sampling

* $\alpha = 0, \beta = 0$
* $N = 20{,}000$

Edge probabilities:

* $\pi_{11} = 0.20$
* $\pi_{12} = 0.05$
* $\pi_{22} = 0.10$

Sampling:

* $p_N \in {0.01, 0.05}$
* Expected samples: 200, 1000

---

### (S.2) Sparse Network & Dense Sampling

* $\beta = \frac{1}{2}, \alpha = 0$
* $N = 100{,}000$

Edge probabilities:

* $\pi_{11} = 0.00506$
* $\pi_{12} = 0.00206$
* $\pi_{22} = 0.00901$

Expected samples: 1000, 5000

---

### (S.3) Sparse Network & Sparse Sampling

* $\beta = \frac{1}{2}, \alpha = \frac{1}{4}$
* $N = 200{,}000$
* $p_N = 0.0002$

Expected samples: ~845

Only **known-parameter case** evaluated

---

### (S.4) Extremely Sparse Model

* Focus: Edge subgraph $H = K_2$
* $\beta = 2, \alpha = 0$
* $N = 1{,}000{,}000$

Edge probabilities:

* $3 \times 10^{-7}, 1 \times 10^{-7}, 2 \times 10^{-7}$

Sampling:

* $p_N \in {0.01, 0.05}$

Quantiles estimated via Monte Carlo (Poisson limit)

---

### (S.5) Dense Model & Extremely Sparse Sampling

* $\beta = 0, \alpha = 1$
* $N = 2{,}000{,}000$
* $p_N = 50/N$

Expected samples: 50

* Ego-centric sampling only

**Results:**

* Empirical coverage: **0.935**
* Average interval width: **0.057**

---

## Key Takeaways

* Prediction intervals depend critically on:

  * Accurate variance estimation
  * Correct limiting distribution (Gaussian vs Poisson)

* Sampling design introduces dependence, complicating inference

* Empirical performance varies significantly across sparsity regimes

---

## References

* Bethe-Hessian clustering: Saade et al. (2014)
* SBM theory and asymptotics: Refer to cited sections in original paper

---

## Appendix

* Variance formulas: Appendix E
* Clustering coefficient variance: Section on clustering coefficient

Additional subgraph cases:

* $K_{1,R-1}$
* $K_R$

---

## Usage

This README is intended for:

* Academic research documentation
* Simulation study replication
* Understanding asymptotic inference in network sampling

---

## Author Notes

Assumes familiarity with:

* Random graph theory
* Stochastic Block Models (SBM)
* Asymptotic statistics

---

**End of README**
