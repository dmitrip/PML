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

## Matlab

Computing the approximate PML distribution for a sample.  Samples must be positive-integer-valued.

    sample = [1 2 1 3 4]'
    empirical_histogram = int_hist(sample) # output: [2 1 1 1]'
    p_PML_approx = PMLdistributionApproximate(empirical_histogram) # output: ones(8,1)./8
    
For the above example, the inferred PML distribution has support size 8, whereas the sample contains only 4 distinct symbols.  If you want to assume a particular support set size, you can supply it as an optinal argument:

    p_PML_approx = PMLdistributionApproximate(empirical_histogram, 4) # output: ones(4,1)./4
 
The approximate PML distribution may not have a finite support set size, in which case we say it has a "continuous" part (Ortlitsky et al. terminology).  In this case, a warning is displayed, the output contains only the discrete (non-continuous) part of the distribution, and the total probability mass of the discrete part of the distribution is less than 1.  For example, the empirical distribution [3 1 1] corresponds to an approximate PML distribution with continuous part of mass 0.4:

    p_PML_approx = PMLdistributionApproximate([3 1 1]) # output: 0.6, a warning is printed

### Estimating functionals of a distribution

The PML distribution can be plugged into a functional of a distribution to obtain the PML plug-in estimate of the functional.  To estimate the entropy and R\'enyi entropy with $\alpha = 2$:

    p_PML_approx = PMLdistributionApproximate([2 1 1 1]) # output: ones(8,1)./8
    H_est = entropyOfDistribution(p_PML_approx) # output: 3 (log base 2)
    Renyi_est = renyiEntropyOfDistribution(p_PML_approx, 2) # output: 3 (log base 2)
    
The functions `estEntroPMLapproximate` and `estRenyiEntroPMLapproximate` wrap the above steps into one function call and compute the PML plug-in estimates given a sample, rather than an empirical histogram.

## Python

Coming soon.

## Julia

Coming soon.
