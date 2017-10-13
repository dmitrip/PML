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

### Estimating a symmetric functional of a distribution 
###### Julia, Matlab, and Python
We first compute the approximate PML distribution "under the hood" and then return the function(al) evaluated on the approximate PML distribution:
```python
F_est = estimate_fun_from_histogram(F, empirical_distribution, [optional] K)
```
where `F` is a function(al) to be estimated and `empirical_distribution` is a collection of non-negative integers.  `K` is an optional argument setting the assumed support set size (must be at least as large as the number of positive entries in `empirical_distribution`).  If `K` is not provided, then we optimize over the support set size.  Zero-valued entries of the empirical distribution are ignored during estimation.  

### Estimating a symmetric functional of multiple distributions 
###### Julia and Python only
If `F` is a function(al) of D distributions -- like L₁ distance for D=2 -- then we need D empirical distributions of the same length to estimate it:
```python
F_est = estimate_fun_from_multiple_histograms(F, [empirical_distribution_1, empirical_distribution_2])
```
This can be used even for a single empirical distribution with D=1 (e.g. estimate entropy), but then you should expect worse performance than using `estimate_fun_from_histogram` from the previous section.  The reason is that for multiple histograms, the PML approximation relies on a heuristic that can be avoided with a special-purpose D=1 implementation.

For now there is no option to set the assumed support set size.

### Computing the approximate PML distribution
###### Julia, Matlab, and Python
```python
p = approximate_PML_from_histogram(empirical_distribution, [optional] K)
```
where `empirical_distribution` is a collection of non-negative integers and `K` is an optional argument setting the assumed support set size (must be at least as large as the number of positive entries in `empirical_distribution`).  If `K` is not provided, then we optimize over the support set size.  Zero-valued entries of the empirical distribution are ignored, so the inferred support size (the length of the output `p`) might be smaller than the length of `empirical_histogram`.

For some inputs, the output `p` has sum less than 1 (for example, if each symbol occurs once, so `empirical_distribution` is a vector of ones).  The missing probability mass is the "continuous part," distributed over infinitely many unobserved symbols, and the output `p` is the "discrete part."

### Computing multiple approximate PML distributions jointly
###### Julia and Python only
Given D empirical distributions of the same length:
```python
p_list = approximate_PML_from_histogram([empirical_distribution_1, empirical_distribution_2, ...])
```
The output `p_list` is a list of the jointly approximated PML distributions.  For now there is no option to set the assumed support set size.

### Examples
###### Python

###### Julia

###### Matlab
