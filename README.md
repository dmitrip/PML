# PML
Profile maximum likelihood approximations and estimation (in Julia, Matlab, and Python)

<p align="center"><img src="https://github.com/dmitrip/PML/blob/master/.github/approx_PML.png" alt="PML approximate distribution" width="50%" alt="PML approximate distribution" /></p>

## Profile maximum likelihood

Suppose we have `n` samples with empirical distribution (histogram) `p̂=(̂p[1], ̂p[2], ...)`.  A relabeling `σ̂p = (p̂[σ[1]], p̂[σ[2]], ...)` permutes the components of `p̂` according to permutation `σ`.  The profile maximum likelihood (PML) distribution maximizes the probability of observing any relabeling of the empirical distribution `p̂`, computed by:
```math
pᴾᴹᴸ = argmaxₚ ∑_σ exp(-n D(σ̂p‖p)) / 𝓕₀!
```
where the sum is over all permutations `σ` of the support set of distribution `p`, `𝓕₀` is the number of symbols seen 0 times empirically, and `D(·‖·)` is the Kullback-Leibler divergence.

The PML distribution can be used as a plug-in estimator for symmetric functionals of a distribution (that is, functionals that are invariant under relabeling, like entropy, Rényi entropy, distance to uniformity, support set size, support set coverage, and others) or functionals of multiple distributions (like L₁ distance, Kullback-Leibler divergence, and others).  [Acharya, Das, Orlitsky, and Suresh 2016](https://arxiv.org/abs/1611.02960) show that the PML approach yields a competitive estimator for symmetric functionals.

The PML distribution is hard to compute, but we can compute it efficiently approximately.  This package implements the approximations presented in [Pavlichin, Jiao, and Weissman 2017]. 


## Usage
Julia, Matlab, and Python implementations share the same interface.  See language-specific examples below.

### Estimating symmetric functionals of distributions

We first compute the approximate PML distribution "under the hood" and then return the function(al) evaluated on the approximate PML distribution.

When the underlying support set size is unknown:
```python
F_est = estimate_fun_from_histogram(F, empirical_distribution)
```
where `F` is a function(al) to be estimated and `empirical_distribution` is a collection of non-negative integers.  Zero-valued entries of the empirical distribution are ignored during estimation.  

When the support set size is assumed to be integer `K` (must be at least as large as the number of positive entries in `empirical_distribution`):
```python
F_est = estimate_fun_from_histogram(F, empirical_distribution, K)
```

### Computing the PML distribution
When the support set size is unknown, then we optimize over it.  Zero-valued entries in `empirical_histogram` are ignored, so the inferred support size (the length of the output `PML_approx`) might be smaller than the length of `empirical_histogram`:
```python
PML_approx = PML_distribution_approximate(empirical_distribution)
```
For some inputs, the output `PML_approx` has sum less than 1 (for example, if each symbol occurs once, so `empirical_distribution` is a vector of ones).  The missing probability mass is the "continuous part," distributed over infinitely many unobserved symbols, and `PML_approx` is the "discrete part."

When the support set size is assumed to be integer `K` (must be at least as large as the number of positive entries in `empirical_distribution`):
```python
PML_approx = PML_distribution_approximate(empirical_distribution, K)
```

