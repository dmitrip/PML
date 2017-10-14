# PML
Approximate profile maximum likelihood estimation.  This package implementats the algorithms in [Pavlichin, Jiao, and Weissman 2017] in Julia, Matlab, and Python.

<p align="left"><img src="https://github.com/dmitrip/PML/blob/master/.github/approx_PML.png" alt="PML approximate distribution" width="50%"/></p>

 * [Profile maximum likelihood overview](#profile-maximum-likelihood-overview)
 * [Performance plots](#performance-plots)
 * [Usage](#usage)
 * [Examples](#examples)
   * [Python](#python)
   * [Julia](#julia)
   * [Matlab](#matlab)
 * [License](#license)

## Profile maximum likelihood overview

Suppose we have `n` samples with empirical distribution (histogram) `p̂=(̂p[1], ̂p[2], ...)`.  A relabeling `σ̂p = (p̂[σ[1]], p̂[σ[2]], ...)` permutes the components of `p̂` according to permutation `σ`.  The profile maximum likelihood (PML) distribution `pᴾᴹᴸ` maximizes the probability of observing any relabeling of the empirical distribution `p̂`, computed by:

<p align="center"><img src="https://github.com/dmitrip/PML/blob/master/.github/eq.png" alt="PML distribution" width="35%"/></p>

where the sum is over all permutations `σ` of the support set of distribution `p`, `𝓕₀` is the number of symbols seen 0 times empirically, and `D(·‖·)` is the Kullback-Leibler divergence.  The support set of `p` is generally not assumed known.  Any relabeling of `pᴾᴹᴸ` maximizes the expression above, so we order `pᴾᴹᴸ` arbitrarily.

The PML distribution can be used in a plug-in estimator `F(pᴾᴹᴸ)` for a symmetric functional `F` of a distribution (that is, a functional that is invariant under relabeling, like entropy, Rényi entropy, distance to uniformity, support set size, support set coverage, and others).  The PML distribution can be generalized to 2 or more distributions approximated jointly, which can be used in a plug-in estimator `F(pᴾᴹᴸ, qᴾᴹᴸ)` for a symmetric functional of multiple distributions (like L₁ distance, Kullback-Leibler divergence, and others).  Acharya, Das, Orlitsky, and Suresh "A Unified Maximum Likelihood Approach for Optimal Distribution Property Estimation" 2016 [link](https://arxiv.org/abs/1611.02960) show that the PML approach yields a competitive estimator for symmetric functionals.

The PML distribution is hard to compute, but we can compute it efficiently approximately.  This package implements the approximations presented in [Pavlichin, Jiao, and Weissman 2017]. 

## Performance plots

Root mean squared error (RMSE) in estimating entropy (unknown support set size) and and L₁ distance from the uniform ditribution (known support set size).  The true support set for the uniform and mixture of 2 uniforms distributions is 10^4, where the mixture of 2 uniforms puts half its mass on the first 20% of symbols, and the other half on the remaining 80% of symbols.  The geometric distribution has mean 10^4.

+ `MLE` is the "naive" empirical distribution plug-in estimator
+ `VV` is Valiant and Valiant, "Estimating the unseen: improved estimators for entropy and other properties" 2017. [link](http://theory.stanford.edu/~valiant/papers/unseenJournal.pdf)
+ `JVHW` is Jiao, Venkat, Han, and Weissman, "Minimax Estimation of Functionals of Discrete Distributions" 2015 [link](https://arxiv.org/abs/1406.6956)
+ `WY` is Wu and Yang, "Minimax rates of entropy estimation on large alphabets via best polynomial approximation" 2014 [link](https://arxiv.org/abs/1407.0381).  This estimator requires the support set size as an input.
+ `approx. PML` is the method implemented here.

<p align="center"><img src="https://github.com/dmitrip/PML/blob/master/.github/approx_PML_performance_github.png" alt="performance of PML approximate distribution plug-in estimator" width="100%"/></p>

RMSE in estimating the L₁ distance between distributions `p` and `q`.  `Zipf(α)` is the distribution `p[i] ~ i^α` up to normalization with finite support 10^4.

<p align="center"><img src="https://github.com/dmitrip/PML/blob/master/.github/approx_PML_performance_L1_distance_github.png" alt="performance of PML approximate distribution plug-in estimator" width="100%"/></p>

## Usage

Julia, Matlab, and Python implementations share the same interface.  See language-specific examples below.

### Estimating a symmetric functional of a distribution 
###### Julia, Matlab, and Python
```python
F_est = estimate_fun_from_histogram(F, empirical_distribution, [optional] K)
```
where `F` is a function(al) to be estimated and `empirical_distribution` is a collection of non-negative integers.  `K` is an optional argument setting the assumed support set size (must be at least as large as the number of positive entries in `empirical_distribution`).  If `K` is not provided, then we optimize over the support set size.  Zero-valued entries of the empirical distribution are ignored during estimation.

We first compute the approximate PML distribution "under the hood" and then return the function(al) evaluated on the approximate PML distribution.  To get the approximate PML distribution, see [below](#computing-the-approximate-pml-distribution).

### Estimating a symmetric functional of multiple distributions 
###### Julia and Python only
If `F` is a function(al) of D distributions -- like L₁ distance for D=2 -- then we need D empirical distributions of the same length to estimate it:
```python
F_est = estimate_fun_from_multiple_histograms(F, [empirical_distribution_1, empirical_distribution_2])
```
This can be used even for a single empirical distribution with D=1 (e.g. estimate entropy), but then you should expect worse performance than using `estimate_fun_from_histogram` from the previous section.  The reason is that for multiple histograms, the PML approximation relies on a heuristic that can be avoided with a special-purpose D=1 implementation.

For now there is no option to set the assumed support set sizes of the D distributions.

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
p_list = approximate_PML_from_multiple_histograms([empirical_distribution_1, empirical_distribution_2, ...])
```
The output `p_list` is a list of the jointly approximated PML distributions.  For now there is no option to set the assumed support set sizes of the D distributions.

### Examples
###### Python
Requires numpy and scipy.  Empirical histogram can be a list or numpy array.
```python
>>> import pml as pml
>>> pml.approximate_PML_from_histogram([2, 1, 1, 1]) # array([ 0.125,  0.125,  0.125,  0.125,  0.125,  0.125,  0.125,  0.125])
```
To assume a particular support set size:
```python
>>> pml.approximate_PML_from_histogram([2, 1, 1, 1], 5) # array([ 0.2,  0.2,  0.2,  0.2,  0.2])
>>> pml.approximate_PML_from_histogram([2, 1, 1, 1], 3) # error, assumed support set size must be at least number of non-zero entries in empirical distribution
```
Some functions of distributions are provided for convenience, others we can define on the fly:
```python
>>> H = pml.entropy_of_distribution # Shannon entropy, log base 2
>>> Renyi = lambda p: pml.renyi_entropy_of_distribution(p, alpha=1.5) # Rényi entropy with α = 1.5, log base 2
>>> support_set_size = lambda p: sum(x > 0 for x in p if x > 0)
>>> L1 = pml.L1_distance # L₁ distance
>>> D_KL = lambda p,q: pml.KL_divergence(p, q, min_ratio=1e-6) # KL divergence with assumed min_x p[x]/q[x] = 1e-6 to avoid infinities
```
Now let's estimate them:
```python
>>> empirical_distribution_1 = [10, 5, 2, 1, 1, 1]
>>> empirical_distribution_2 = [2, 1, 1, 0, 1, 1] # must be same length as empirical_distribution_1 to estimate L₁ distance
>>> pml.estimate_fun_from_histogram(H, empirical_distribution_1) # 2.25
>>> pml.estimate_fun_from_histogram(Renyi, empirical_distribution_1) # 1.8716...
>>> pml.estimate_fun_from_histogram(support_set_size, empirical_distribution_1) # 10
>>> pml.estimate_fun_from_multiple_histograms(L1, [empirical_distribution_1, empirical_distribution_2]) # 
>>> pml.estimate_fun_from_multiple_histograms(D_KL, [empirical_distribution_1, empirical_distribution_2]) # 
```
###### Julia

###### Matlab

### License
This project is licensed under the [Apache License 2.0](LICENSE.md).
