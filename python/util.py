import numpy as np

def draw_histogram_from_multinomial(p, n):
    """
    draw histogram on n samples from multinomial distribution specified by weight vector p (need not be normalized)

    Args:
        p (1d array): weights
        n (int): sample size

    Returns:
        q (1d array int): sample histogram

    """

    p = p/np.sum(p)    
    q = np.zeros(p.shape, dtype=int) # current type
    z = np.sum(q)
    num_trials = 0
    alpha = 0.99
    while z < n:
        q_ = np.random.poisson(alpha*(n-z)*p)
        z_ = np.sum(q_)
        if z + z_ <= n: # accept
            q += q_
            z += z_
        num_trials += 1
    return q

