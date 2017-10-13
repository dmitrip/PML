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

## Estimating symmetric functionals of distributions

To estimate a functional of a probability distribution we first compute the approximate PML distribution "under the hood" and then return the functional evaluated on the approximate PML distribution.

When the support set size is unknown:
```python
F_est = estimate_fun_from_histogram(F, empirical_distribution)
```
where `F` is a function(al) to be estimated and `empirical_distribution` is a collection of non-negative integers.  Zero-valued entries of the empirical distribution are ignored during estimation.  

To assume a particular value for the support set size:
```python
F_est = estimate_fun_from_histogram(F, empirical_distribution, K)
```
where `K` is the integer-valued assumed support set size, must be at least as large as the number of positive entries in `empirical_distribution`.

## Computing the PML distribution

