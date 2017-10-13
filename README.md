# PML
Profile maximum likelihood approximations

<a name="logo"/>
<div align="center">
<img src="https://github.com/dmitrip/PML/blob/master/.github/PML_approximation.png" alt="PML approximate distribution"></img>
</a>
</div>

## Profile maximum likelihood

For samples $x_1^n = (x_1,\ldots,x_n)$ with empirical distribution $\hat{p}$, the profile maximum likelihood (PML) distribution $p^*_\text{PML}$ maximizes the probability of observing any relabeling of $\hat{p}$:

$$p^*_\text{PML} = \arg \max_p \sum_{\sigma \in S_\mathcal{X}} \mathbb{P}_p(\sigma \hat{p})$$

While $p^*$ is hard to compute (the optimization involves maximizing a matrix permanent), we can compute it efficiently approximately.  This package implements approximate PML and exact PML in small cases.

The PML distribution can be used as a plug-in estimator for symmetric functionals of distributions (that is, functionals that are invariant under relabeling of the support set, like entropy, Rényi entropy, support set size).

## Usage
The Julia, Matlab, and Python implementations share the same interface.  See examples below.

### Estimating symmetric functionals of distributions

To estimate a functional of a probability distribution we first compute the approximate PML distribution "under the hood" and then return the functional evaluated on the approximate PML distribution.

When the support set size is unknown:
```python
F_est = estimate_fun_from_histogram(F, empirical_distribution)
```
where `F` is a function(al) to be estimated and `empirical_distribution` is a collection of non-negative integers.  Zero-valued entries of the empirical distribution are ignored during estimation.  

When the support set size is assumed to be integer `K` (must be at least as large as the number of positive entries in `empirical_distribution`):
```python
F_est = estimate_fun_from_histogram(F, empirical_distribution, K)
```

### Computing the PML distribution
If the underlying support set size is unknown, then we optimize over it.  Zero-valued entries in the `empirical_histogram` are ignored, so the inferred support might be smaller than the length of `empirical_histogram`:
```python
PML_approx = PML_distribution_approximate(empirical_distribution)
```
It's possible that the output `PML_approx` has sum less than 1 (for example, if each symbol occurs once, so `empirical_distribution` is a vector of ones).  The missing probability mass is the "continuous part" [Orlitsky, Santhanam, Viswanathan, Zhang, 2004], distributed over infinitely many symbols, and `PML_approx` is the "discrete part."

When the support set size is assumed to be integer `K` (must be at least as large as the number of positive entries in `empirical_distribution`):
```python
PML_approx = PML_distribution_approximate(empirical_distribution, K)
```

