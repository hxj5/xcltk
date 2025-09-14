# hmm.py
# - copied from Github repo `cfDNA2karyotype` (Author Yuanhua Huang).


import numpy as np
from scipy.special import logsumexp


def logdotexp(A_log, B_log):
    """
    Extension of logsumexp to logdotexp
    Note, this is not memory efficient and will have peak size: (n, k, m) 
    for (n, k) * (k, m)
    
    Examples
    --------
    >>> A = np.arange(1, 7).reshape(3, 2)
    >>> B = np.arange(1, 9).reshape(2, 4)
    >>> print(np.log(A @ B))
    >>> print(np.log(np.dot(A, B)))
    >>> print(logdotexp(np.log(A), np.log(B)))
    """
    if len(A_log.shape) == 1:
        A_log = A_log.reshape(1, -1)
    if len(B_log.shape) == 1:
        B_log = B_log.reshape(-1, 1)
        
    AB_log = logsumexp(
        np.expand_dims(A_log, 2) + np.expand_dims(B_log, 0),
        axis=1
    )
    
    return AB_log



def HMM(emm_p_log, pi=None, A=None, diag_prob=0.99):
    """
    Hidden Markov Model
    -------------------
    Forward-Backward algorithm for calculating the 
    posterior of the state assignment via HMM
    """
    n_index = emm_p_log.shape[0]
    n_state = emm_p_log.shape[1]

    # prior distribution
    if pi is None:
        pi = np.ones(n_state) / n_state

    # transition matrix
    if A is None:
        diag_prob_adj = (diag_prob - 1 / n_state) * n_state / (n_state - 1)
        A = np.eye(n_state) * diag_prob_adj + (1 - diag_prob_adj) / n_state
    A_log = np.log(A)
    
    # Forward-Backward algorithm for updating posterior
    fw_p_log = emm_p_log.copy()
    bw_p_log = emm_p_log.copy()

    # Forward part of the algorithm
    fw_p_log[0, :] += np.log(pi)
    for i in range(1, n_index):
        fw_p_log[i, :] += logdotexp(fw_p_log[i - 1, :], A_log).reshape(-1)

    # Backward part of the algorithm
    bw_p_log[-1, :] = 0
    for i in range(n_index - 1, 0, -1):
        bw_p_log[i - 1, :] = logdotexp(
            bw_p_log[i, :] + emm_p_log[i, :], A_log).reshape(-1)

    # Update posterior of state assignment
    z_post_log  = fw_p_log + bw_p_log
    z_post_log -= logsumexp(z_post_log, axis=1, keepdims=True)
    z_post = np.exp(z_post_log)
    z_post = z_post / np.sum(z_post, axis=1, keepdims=True)
    
    return z_post
