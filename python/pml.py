import numpy as np
np.seterr(invalid='ignore') 
np.seterr(divide='ignore') # for case 0/0
from scipy import special # for gammaln
import warnings

def estimate_entropy_PML_approximate(samp,K=None):
    # estimates entropy of distribution given sample samp
    # optionally given assumed true support set size K
    #
    # approximates pattern maximum likelihood (PML) distribution p
    # where p \approx \arg \max_p' P(unordered histogram)
    # Approximation sums over all permutations that mix within blocks of
    # constant p
    #
    # Uses approximate PML distribution as plug-in for entropy functional
    #
    # created by Dmitri S. Pavlichin on May 31, 2017
    #
    # Args:
    #     * samp - (vector) sample, vector of non-negative integers
    # Optional args:
    #     * K - (integer) assumed support set size.  
    #         If K is not provided, then we attempt to estimate the support set size
    #
    # Returns:
    #     * H_est - (double) estimate for entropy, log base 2
    #     * K_est - (integer) estimate for support set size
    #
    # Example:

    # get empirical histogram
    hist_vec = int_hist(samp)

    if K is not None:
        K_est = K
        p_approx = PMLdistributionApproximate(hist_vec,K)
        not_have_valid_K_est = False
    else:
        # estimate support set size
        (p_approx,F0_est,_,not_have_valid_K_est) = PML_distribution_approximate(hist_vec)
        K_est = np.sum(hist_vec > 0) + F0_est

    if not_have_valid_K_est:
        warnings.warn('do not have valid estimate for support set size')

    # plug-in to entropy functional
    H_est = entropy_of_distribution(p_approx, 2)

    return H_est

def estimate_Renyi_entropy_PML_approximate(samp,alpha,K=None):
    # estimates Renyi entropy with parameter alpha of distribution given sample samp
    # optionally given assumed true support set size K
    #
    # approximates pattern maximum likelihood (PML) distribution p
    # where p \approx \arg \max_p' P(unordered histogram)
    # Approximation sums over all permutations that mix within blocks of
    # constant p
    #
    # Uses approximate PML distribution as plug-in for Renyi entropy functional
    #
    # created by Dmitri S. Pavlichin on May 31, 2017
    #
    # Args:
    #     * samp - (vector) sample, vector of non-negative integers
    #     * alpha - (double) Renyi entropy parameter, alpha >= 0, alpha != 1
    # Optional args:
    #     * K - (integer) assumed support set size.  
    #         If K is not provided, then we attempt to estimate the support set size
    #
    # Returns:
    #     * H_est - (double) estimate for Renyi entropy, log base 2
    #     * K_est - (integer) estimate for support set size
    #
    # Example:

    # get empirical histogram
    hist_vec = int_hist(samp)

    if K is not None:
        K_est = K
        p_approx = PMLdistributionApproximate(hist_vec,K)
        not_have_valid_K_est = False
    else:
        # estimate support set size
        (p_approx,F0_est,_,not_have_valid_K_est) = PML_distribution_approximate(hist_vec)
        K_est = np.sum(hist_vec > 0) + F0_est

    if not_have_valid_K_est:
        warnings.warn('do not have valid estimate for support set size')

    # plug-in to Renyi entropy functional
    RenyiH_est = renyi_entropy_of_distribution(p_approx, alpha, 2)

    return RenyiH_est


def entropy_of_distribution(vec,base=2.0):
    # computes Shannon entropy of vector vec, default base 2
    Z = np.sum(vec)
    H = sum(-p*np.log(p) for p in vec if p > 0.0)
    return (H/Z + np.log(Z))/np.log(base)

def renyi_entropy_of_distribution(vec, alpha, base=2.0):
    # computes Renyi entropy of vector vec, parameter alpha, default base 2
    Z = np.sum(vec)
    return (np.log(np.sum(np.power(p/Z,alpha) for p in vec))/(1-alpha)) / np.log(base)

def int_hist(sample):
    """get integer histogram of all integer values in range(np.max(sample))
    
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
    C[i,j] = sum(x[i:j+1])

    """    
    K = np.size(x)
    C = np.zeros((K,K))
    for i in range(K):
        C[i,i] = x[i]
        for j in range(i+1,K):
            C[i,j] = C[i,j-1] + x[j]
    
    return C

def PML_distribution_approximate(p_empirical, K=None):
    # approximates pattern maximum likelihood (PML) distribution p_approx
    # where p_approx \approx \arg \max_p P_p(unordered p_empirical)
    # Approximation sums over all permutations that mix within blocks of
    # constant p_approx
    #
    # created by Dmitri S. Pavlichin on May 14, 2017
    #
    # Args:
    #     * p_empirical - (integer-valued vector) empirical histogram, entries sum to sample size
    # Optional args:
    #     * K - (integer) assumed support set size, must have K >= sum(p > 0).  
    #         If K is not provided, then we attempt to estimate the support set size
    # #     * K_upper_bound - (integer) assumed maximum value for support set
    #         size, use when estimating support set size.  If K_upper_bound is
    #         not provided, then we set K_upper_bound = sum(p == 1)^2
    #
    # Returns:
    #     * p_approx - (column vector) approximate PML distribution, sorted in
    #         descending order
    #     * F0_est - (integer) estimate for number of symbols seen 0 times in
    #         p_empirical.  If K is provided as an argument, then F0_est = max(K - sum(p > 0), 0)
    #     * V_approx - (double) approximate value of log(P_p(unordered p_empirical))
    #         achieved by p = p_approx
    #     * reached_max_F0 - (boolean) true iff we did not find an upper bound
    #        when estimating support set size, hit upper bound sum(p > 0) + sum(p == 1)^2
    #
    # Example:
    #     Optimize support set size:
    #         >> p_approx = PMLdistributionApproximate([9,3,2,1,1])
    #         yields p_approx == [9/16 7/80 7/80 7/80 7/80 7/80]'
    #     OK to have zeros in input histogram, zeros are ignored:
    #         >> p_approx = PMLdistributionApproximate([9 3 2 1 1 0 0 0 0 0 0]')
    #         yields p_approx == [9/16 7/80 7/80 7/80 7/80 7/80]'
    #     Assumed support set size:
    #         >> p_approx = PMLdistributionApproximate([9 3 2 1 1]',5)
    #         yields p_approx == [9/16 7/64 77/6/64 7/64 7/64]'

    def f0(V, counts, multiplicities, F0, n, num_bins, B):
        Vmax = special.gammaln(num_bins+F0+1) + n*np.log(1/(num_bins+F0))
        imax = 0
        for i in range(B):
            c = np.sum(counts[i+1:])
            m = np.sum(multiplicities[i+1:])

            if c > 0:
                proposal_reward = V[i] + special.gammaln(m+F0+1) + c*np.log((c/n)/(m+F0))
            else:
                proposal_reward = V[i] + special.gammaln(m+F0+1)

            if proposal_reward > Vmax:
                Vmax = proposal_reward
                imax = i+1

        return (Vmax,imax)

    p_empirical = np.array(p_empirical)
    
    fing = int_hist(p_empirical[p_empirical>0])
    bin_values = np.flipud(np.where(fing>0)[0]) # sort in descending order
    multiplicities = fing[bin_values]
    num_bins = np.sum(multiplicities)

    B = bin_values.size # multibins
    n = np.dot(bin_values, multiplicities)
    n_1 = np.sum(p_empirical==1) # number of symbols observed once

    # early detect case n_1 == n?

    # get counts and multiplicities cumsum matrix
    # n[i,j] = sum_{k = i to j} n[k]
    counts = bin_values*multiplicities
    counts_mat = cumsum_mat(counts)
    mults_mat = cumsum_mat(multiplicities)
    prob_mat = counts_mat/(n*mults_mat) # probability within each block

    log_reward_mat = special.gammaln(mults_mat+1) + counts_mat*np.log(prob_mat)
    log_reward_mat[np.isnan(log_reward_mat)] = 0.0

    # dynamic programming
    # V[i] = best log reward for first i bins
    V = np.zeros(B) # % V[-1] = 0, values matrix
    backpointers = np.zeros(B, dtype=int) # backpointers to get best log reward
    for i in range(B):
        V[i] = log_reward_mat[0,i]
        for j in range(i):
            proposal_reward = V[j] + log_reward_mat[j+1,i]
            if proposal_reward > V[i]:
                V[i] = proposal_reward
                backpointers[i] = j+1

    ## optimize alphabet size
    # get upper bound on F0
    reached_max_F0 = False
    if K is None:
        F0 = 0 # extra 0 entries
        step_size = 1
        F0_max = num_bins**2
        V0 = f0(V, counts, multiplicities, 0, n, num_bins, B)[0] - special.gammaln(0+1)

        done = False
        while not done:
            F0 += step_size

            V0_prev = V0
            V0 = f0(V, counts, multiplicities, F0, n, num_bins, B)[0] - special.gammaln(F0+1)

            if V0_prev >= V0:
                done = True
            else:
                step_size *= 2
                if F0 > F0_max:
                    reached_max_F0 = True
                    done = True
        F0_max = F0

        # get optimum value of F0
        if not reached_max_F0:
            F0_lower = 0
            F0_upper = F0_max
            while F0_lower != F0_upper:
                F1 = (F0_upper-F0_lower)//3+F0_lower
                F2 = 2*(F0_upper-F0_lower)//3+F0_lower+1
                V1 = f0(V, counts, multiplicities, F1, n, num_bins, B)[0] - special.gammaln(F1+1)
                V2 = f0(V, counts, multiplicities, F2, n, num_bins, B)[0] - special.gammaln(F2+1)
                if V1 - V2 <= -1e-10:
                    F0_lower = F1+1
                else:
                    F0_upper = F2-1
            F0 = F0_lower
        else:
            F0 = F0_max

        # test to see if have continuous part
        if bin_values[-1] == 1:
            V_non_cont = f0(V, counts, multiplicities, F0, n, num_bins, B)[0] - special.gammaln(F0+1)
            if bin_values.size > 1:
                V_cont = V[-2] + n_1*np.log(n_1/n)
            else:
                V_cont = 0.0
            if V_cont > V_non_cont:
                warnings.warn('optimal distribution has continuous part, F0 = inf')
                F0 = np.inf
    else:
        F0 = K - num_bins
        if F0 < 0:
            F0 = 0
    F0_est = F0

    if not np.isinf(F0):
        (Vmax,imax) = f0(V, counts, multiplicities, F0, n, num_bins, B)
        backpointers = np.append(backpointers, imax)
        bin_values = np.append(bin_values, 0)
        multiplicities = np.append(multiplicities, F0)
        B += 1
        num_bins += F0
        counts = bin_values*multiplicities
        counts_mat = cumsum_mat(counts)
        mults_mat = cumsum_mat(multiplicities)
        prob_mat = counts_mat/(n*mults_mat) # probability within each block
    else:
        Vmax = V_cont
        backpointers = backpointers[:-1]
        bin_values = bin_values[:-1]
        multiplicities = multiplicities[:-1]
        B -= 1
        num_bins -= n_1
        counts = bin_values*multiplicities
        counts_mat = cumsum_mat(counts)
        mults_mat = cumsum_mat(multiplicities)
        prob_mat = counts_mat/(n*mults_mat) # probability within each block

    ## get assignment
    p_approx = np.zeros(num_bins)

    current_multibin = B-1
    ix = num_bins-1

    while ix >= 0:
        t_prev = current_multibin
        t = backpointers[current_multibin]
        while current_multibin >= t:
            for m in range(multiplicities[current_multibin]):
                p_approx[ix] = prob_mat[t,t_prev]
                ix -= 1
            current_multibin -= 1

    ##

    V_approx = Vmax + special.gammaln(n+1) \
        - np.sum(multiplicities*special.gammaln(bin_values+1)) \
        - np.sum(special.gammaln(multiplicities+1))
    
    return (p_approx, F0_est, V_approx, reached_max_F0)

