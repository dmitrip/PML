# created by Dmitri S. Pavlichin in October 2017

import numpy as np
np.seterr(invalid='ignore') 
#np.seterr(divide='ignore') # for case 0/0
from scipy import special # for gammaln
#import warnings
import collections # for Counter

def counts(sample):
    """get histogram of sample, ordered arbitrarily
    
    Args:
        sample - list
    
    Returns:
        numpy array h, where h[i] = |{t: sample[t] == i}|
        
    Examples:
        >>> counts([0,0,1,'cat',1,1])
        [2, 3, 1]
    
    """
    
    _, h = np.unique(sample, return_counts=True)
    return h

def int_hist(sample):
    """get integer histogram of all integer values in range(np.max(sample)+1)
    
    Args:
        sample - integer-valued list
    
    Returns:
        h - h[i] = |{t: sample[t] == i}|
        
    Examples:
        >>> int_hist([0,0,1,5,1,1])
        array([2, 3, 0, 0, 0, 1])
    
    """
    sample_max = np.max(sample)
    return np.histogram(sample, bins=sample_max+1, range=(0,sample_max))[0]

def cumsum_mat(x):
    """computes cumultative sum matrix
    x is numpy array
    C[i,j] = sum(x[i:j+1])

    """    
    K = np.size(x)
    C = np.zeros((K,K), dtype=x.dtype)
    for i in range(K):
        C[i,i] = x[i]
        for j in range(i+1,K):
            C[i,j] = C[i,j-1] + x[j]
    
    return C

def maximize_on_interval_int(fun, x_min, x_max, tol=1e-14):
    """maximizes semiconcave function fun, integer arguments, on interval [x_min, x_max]

    created by Dmitri S. Pavlichin on May 14, 2017

    Args:
        * fun - function to maximize
        * x_min - interval min
        * x_max - interval max

    Optional:
        * tol (float) - tolerance optimize

    Assumptions:
        * interval [x_min, x_max] includes at least one integer

    Returns:
        (x, val)
            * x (int) - maximizing value of x
            * val - fun(x)

    """

    # check assumptions
    assert np.ceil(x_min) <= np.floor(x_max), '[x_min, x_max] interval must include at least one integer'

    while np.floor(x_max) > np.floor(x_min):
        x1 = (x_max - x_min)/3 + x_min
        x2 = 2*(x_max - x_min)/3 + x_min
        V1 = fun(x1)
        V2 = fun(x2)
        if V1 - V2 <= -np.abs(tol):
            x_min = x1
        else:
            x_max = x2

    x = x_min
    if fun(np.floor(x)) >= fun(np.ceil(x)):
        x = np.floor(x_min)
    else:
        x = np.ceil(x_min)

    return (x, fun(x))

def harmonic_number(n):
    # harmonic number H_n
    return np.sum(1/i for i in range(1,n+1))

def ML_unseen_symbols_uniform(T, n):
    """maximum likelihood estimate for number of unseen symbols given empirical
    support set size T, n samples

    created by Dmitri S. Pavlichin on October 2, 2017

    Args:
        * T (int) - empirical support set size
        * n (int) - sample size

    Returns:
        * F0 (int) - maximum likelihood estimator for number of unseen symbols

    """

    # check for F0 == 0
    if n >= T*harmonic_number(T):
        F0 = 0
    else:
        # check for F0 == inf
        if T >= n:
            F0 = np.inf
        else:
            fun = lambda f0: special.gammaln(f0+T+1)-special.gammaln(f0+1)-n*np.log(f0+T)
            F0 = maximize_on_interval_int(fun, 0, np.ceil((T*T-n)/(n-T))+1)[0]

    return F0

def approximate_PML_from_histogram(p_empirical, K=None):
    """
    approximates pattern maximum likelihood (PML) distribution p_approx
    where p_approx \approx \arg \max_p P_p(unordered p_empirical)
    Approximation sums over all permutations that mix within blocks of
    constant p_approx
    
    created by Dmitri S. Pavlichin on May 14, 2017
    
    Args:
        * p_empirical - (integer-valued vector) empirical histogram, entries sum to sample size
    Optional args:
        * K - (integer) assumed support set size, must have K >= sum(p > 0).  
            If K is not provided, then we attempt to estimate the support set size
    
    Returns:
        * p_approx - (column vector) approximate PML distribution, sorted in
            descending order
        * F0_est - (integer) estimate for number of symbols seen 0 times in
            p_empirical.  If K is provided as an argument, then F0_est = max(K - sum(p > 0), 0)
        * V_approx - (double) approximate value of log(P_p(unordered p_empirical))
            achieved by p = p_approx
    
    Example:
        Optimize support set size:
            >>> p_approx = approximate_PML_from_histogram([9,3,2,1,1])
            returns p_approx == [9/16 7/80 7/80 7/80 7/80 7/80]'
        OK to have zeros in input histogram, zeros are ignored:
            >>> p_approx = approximate_PML_from_histogram([9 3 2 1 1 0 0 0 0 0 0]')
            returns p_approx == [9/16 7/80 7/80 7/80 7/80 7/80]'
        Assumed support set size:
            >>> p_approx = approximate_PML_from_histogram([9 3 2 1 1]',5)
            returns p_approx == [9/16 7/64 77/6/64 7/64 7/64]'

    """

    p_empirical = np.array(p_empirical)
    p_empirical = p_empirical[p_empirical != 0] # remove 0-counts
    
    fing = int_hist(p_empirical)
    fing_supp = np.where(fing>0)[0] # sorted in ascending order
    fing_values = fing[fing_supp]

    B = fing_supp.size # multibins, F_+
    n = np.dot(fing_supp, fing_values) # bin counts times fing_values, sum(bm) = sum(p)

    # early detect case all symbols observed once?

    # get counts and multiplicities cumsum matrix
    # n[i,j] = sum_{k = i to j} n[k]
    counts = fing_supp*fing_values
    counts_mat = cumsum_mat(counts)
    mults_mat = cumsum_mat(fing_values)
    prob_mat = counts_mat/(n*mults_mat) # probability within each block

    # C[i,j] = log reward assign multibins from i to j to same prob
    log_reward_mat = special.gammaln(mults_mat+1) + counts_mat*np.log(prob_mat)
    log_reward_mat[np.isnan(log_reward_mat)] = 0.0 # handles 0 log 0 = 0 case, sets below diagonal entries to 0

    # dynamic programming
    # V[i] = best log reward for first i bins
    V = np.zeros(B) # % V[-1] = 0, values matrix
    backpointers = np.zeros(B, dtype=int)-1 # backpointers to get best log reward
    for i in range(B-1,-1,-1):
        V[i] = log_reward_mat[i,B-1]
        backpointers[i] = B-1
        for j in range(i,B-1):
            proposal_reward = log_reward_mat[i,j] + V[j+1]
            if proposal_reward > V[i]:
                V[i] = proposal_reward
                backpointers[i] = j

    ## merge 0s with rest of symbols
    V0_opt = -np.inf # optimal value
    F0_opt = -1
    i0_opt = -1 # index of optimal fingerprint bin up to which to merge 0s
    for i in range(B):
        T = np.sum(fing_values[:(i+1)])
        N = np.sum(fing_supp[:(i+1)]*fing_values[:(i+1)])
        if K is not None: # assumed alphabet size
            F0 = K - len(p_empirical)
        else:
            # optimize alphabet size
            # estimate number of unseen symbols for uniform distribution, N samples, empirical support set size T
            F0 = ML_unseen_symbols_uniform(T, N)
        # get value
        V0 = N*np.log(N/n)
        if i < B-1:
            V0 += V[i+1] # V[B] = 0
        if not np.isinf(F0): # \lim_{F0 -> \infty} f(F0,T,n=T) = 0
            V0 += (special.gammaln(F0+T+1) - special.gammaln(F0+1) - N*np.log(F0+T))
        # check if best so far
        if V0 > V0_opt:
            V0_opt = V0
            F0_opt = F0
            i0_opt = i
    F0 = F0_opt
    V_approx = V0_opt

    ## get assignment
    p_approx = np.zeros_like(p_empirical, dtype=float)
    if F0 > 0:
        backpointers[0] = i0_opt # merge 1s with same symbols as 0s, since guaranteed to merge 0s at least with 1s 
    i = 0
    while i <= B-1:
        j = backpointers[i]
        if i == 0:
            x1 = 0
        else:
            x1 = mults_mat[0,i-1]
        x2 = mults_mat[0,j]
        p_approx[x1:x2] = prob_mat[i,j]
        i = j + 1
    continuous_part = 0 # probability mass of continuous part
    # add 0s in front
    if F0 > 0:
        x = mults_mat[0,i0_opt] # number of seen symbols in same level set as unseen symbols
        if not np.isinf(F0):
            p_approx = np.insert(p_approx, 0, np.zeros(int(F0)))
            p_approx[0:(x+int(F0))] = counts_mat[0,i0_opt]/(n*(x+F0))
        else:
            continuous_part = counts_mat[0,i0_opt]/n
            p_approx = p_approx[x:]

    ## get log probability of fingerprint
    V_approx = V_approx + special.gammaln(n+1) - np.sum(fing_values*special.gammaln(fing_supp+1)) - np.sum(special.gammaln(fing_values+1))

    ## sort PML approximation in descending order
    p_approx = np.flipud(p_approx)

    have_continuous_part = np.isinf(F0)

    #return (p_approx, F0, V_approx, have_continuous_part)
    return p_approx

def entropy_of_distribution(vec,base=2.0):
    # computes Shannon entropy of vector vec, default base 2
    Z = np.sum(vec)
    H = sum(-p*np.log(p) for p in vec if p > 0.0)
    return (H/Z + np.log(Z))/np.log(base)

def renyi_entropy_of_distribution(vec, alpha, base=2.0):
    # computes Renyi entropy of vector vec, parameter alpha, default base 2
    Z = np.sum(vec)
    return (np.log(np.sum(np.power(p/Z,alpha) for p in vec))/(1-alpha)) / np.log(base)

def estimate_fun_from_histogram(F, p_empirical, K=None):
    # estimate function F of distribution
    p_PML_approx = approximate_PML_from_histogram(p_empirical, K)
    return F(p_PML_approx)


