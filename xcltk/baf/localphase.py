# localphase.py - local phasing
# further check the SNP-level reference phasing states after running Eagle2.


import numpy as np
import scipy as sp
from logging import info
from logging import warning as warn
from scipy.special import logsumexp
from ..utils.hmm import HMM



def snp_local_phasing(
    AD, DP,
    min_iter = 5, max_iter = 50,
    cw_min_expr_snps = 1,
    cw_low_baf = 0.45, cw_up_baf = 0.55,
    positions = None,
    smooth_hmm = False,
    kws_hmm = None,
    smooth_gk = True,
    kws_gk = None,
    kws_local_phasing = None,
    verbose = False
):
    """Phase SNPs within small blocks by assuming the allelic ratio is the same
    for all SNPs and adjacent SNPs are more likely to be flipped together if
    any one is to be flipped.

    To alleviate the effect of noisy non-CNA cells on phasing, an iterative
    strategy is used to remove cells whose BAF (aggregated from all SNPs) is
    near 0.5.

    Parameters
    ----------
    AD : np.array like
        The cell x snp AD matrix.
    DP : np.array like
        The cell x snp DP matrix.
    cw_min_expr_snps : int
        Minimum number of expressed SNPs in one cell. 
        Used for cell filtering in each iteration.
    cw_low_baf : float
        The cell-wise lower bound of aggregated BAF.
        Cells between `cw_low_baf` and `cw_up_baf` will be filtered in each
        iteration.
    cw_up_baf : float
        The cell-wise up bound of aggregated BAF.
    positions : array-like
        The SNP positions.
    smooth_hmm : bool
        Whether to use HMM for haplotype assignment smoothing.
    kws_hmm : dict or None
        Other parameters passed to `HMM()`.
    smooth_gk : bool
        Whether to use gaussian kernel for haplotype assignment smoothing.
    kws_gk : dict or None
        Other parameters passed to `snp_smooth_gk()`.
    kws_local_phasing : dict or None
        Other parameters passed to `Local_Phasing()`.

    Returns
    -------
    None or phasing results.
        None if error or no sufficient cells for phasing.
    """
    if kws_local_phasing is None:
        kws_local_phasing = dict()
    if kws_hmm is None:
        kws_hmm = dict()
        
    idx = DP.sum(axis = 1) > 0
    AD, DP = AD[idx, :], DP[idx, :]
    if AD.shape[0] <= 0:
        if verbose:
            warn("%d - all cells are filtered (low DP)!")
        return None
    
    n_expr_snps = (DP > 0).sum(axis = 1)
    idx = n_expr_snps >= cw_min_expr_snps
    AD, DP = AD[idx, :], DP[idx, :]
    if AD.shape[0] <= 0:
        if verbose:
            warn("%d - all cells are filtered (few expressed SNPs)!")
        return None
        
    flip_final = None
    for i in range(max_iter):
        baf = AD.sum(axis = 1) / DP.sum(axis = 1)
        idx = np.logical_or(baf < cw_low_baf, baf > cw_up_baf)
        AD, DP = AD[idx, :], DP[idx, :]
        if AD.shape[0] <= 0:
            if verbose:
                warn("%d - all cells are filtered (BAF near 0.5)!")
            return None

        ad_sum, ad_sum1, dp_sum, Z, thetas, logLik_new = Local_Phasing(
            AD = AD.T,
            DP = DP.T,
            smooth_hmm = smooth_hmm,
            kws_hmm = kws_hmm,
            positions = positions,
            smooth_gk = smooth_gk,
            kws_gk = kws_gk,
            verbose = verbose,
            **kws_local_phasing
        )
        
        flip = np.array(Z[:, 1] >= Z[:, 0])
        if i == 0:
            flip_final = flip
        else:
            flip_final = np.logical_xor(flip_final, flip)
        if verbose:
            info("%d - %s %f %d" % (i, str(AD.shape), logLik_new, np.sum(flip)))
            info("    flip_final: %s" % str(flip_final))
            info("    flip: %s" % str(flip))
            info("    Z: %s" % str(Z))
        
        if i > 0:
            if i >= min_iter and (np.all(flip) or np.all(np.logical_not(flip))):
                if verbose:
                    info("%d - flip state converged." % i)
                return flip_final, ad_sum, ad_sum1, dp_sum, Z, thetas, logLik_new

        AD = AD * (1 - flip.T) + (DP - AD) * flip.T

    return flip_final, ad_sum, ad_sum1, dp_sum, Z, thetas, logLik_new


                             
# This function was copied from XClone repo.
# Modified by XH.
def Local_Phasing(
    AD, DP, 
    min_iter = 10, max_iter = 1000, epsilon_conv = 1e-3,
    init_mode = 'warm', 
    positions = None,
    smooth_hmm = False,
    kws_hmm = None,
    smooth_gk = False,
    kws_gk = None,
    verbose = False
):
    """
    Phase the small blocks into a medium sized bin by assuming the allelic
    ratio is the same for all blocks. This is equavilent to a binary clustering,
    whose likelihood can be maximised by an EM alogrithm

    # Note that AD DP should be in in Compressed Sparse Column format.
    # or other sparse format.
    
    Parameters
    ----------
    AD : sparse matrix of integers
        Read counts for ALT allele in N blocks and M cells.
    DP : sparse matrix of integers
        Read counts for REF allele in N blocks and M cells.
    
    Returns
    -------
    ad_sum, dp_sum, Z, thetas, logLik
    """
    def __safe_loglik(AD, DP, thetas, eps = 1e-6):
        thetas = thetas.copy()
        thetas[thetas <= 0] = eps
        thetas[thetas >= 1] = 1 - eps
        return AD @ np.log(thetas) + BD @ np.log(1 - thetas)

    if smooth_hmm and kws_hmm is None:
        kws_hmm = dict()
    if smooth_gk:
        if positions is None:
            raise ValueError("[Local_Phasing] coordinates should be specified!")
        if kws_gk is None:
            kws_gk = dict()

    eps = 1e-6        # Extreamly small value, e.g., 1e-20, will produce NaN in `Z` and `thetas`.
    N, M = AD.shape
    BD = DP - AD
    
    ## Initialization matters!!!
    if init_mode == 'warm':        # Option 1: warm initialization (usually good?)
        Z = np.zeros((N, 2))
        Z[:, 0] = (AD.sum(1) / DP.sum(1)).reshape(-1)
        Z[:, 1] = 1 - Z[:, 0]
    elif init_mode == 'current':   # Option 2: initializing with no flipping (can return poor local optimal)
        Z = np.zeros((N, 2))
        Z[:, 0] = 1       # high prob means no flipping
    else:             # Option 3: random initialization (may need large number of trials)
        Z = np.random.rand(N, 2)
        Z[:, 1] = 1 - Z[:, 0]
    
    # allele ratio parameters
    thetas = np.array((AD.T @ Z + BD.T @ (1 - Z)) / (DP.T.sum(1, keepdims = True)))
    
    # likelihood
    _logLik_mat = __safe_loglik(AD, DP, thetas, eps)
    _logLik_new = np.sum(logsumexp(_logLik_mat, axis = 1))
    
    for it in range(max_iter):
        _logLik_old = _logLik_new + 0.0
        
        # E step: calculate the expecation
        Z = normalize(np.exp(loglik_amplify(np.array(_logLik_mat))))
        if smooth_hmm:
            Z = HMM(Z, **kws_hmm)
        if smooth_gk:
            Z[:, 0] = snp_smooth_gk(Z[:, 0], positions, **kws_gk)
            Z[:, 1] = 1 - Z[:, 0]
        
        # M step: maximise the likihood over thetas
        thetas = np.array((AD.T @ Z + BD.T @ (1 - Z)) / (DP.T.sum(1, keepdims = True)))
        
        # Likelihood
        _logLik_mat = __safe_loglik(AD, DP, thetas, eps)
        _logLik_new = np.sum(logsumexp(_logLik_mat, axis = 1))
        
        if verbose:
            info("%d %f %f" % (it, _logLik_old, _logLik_new))

        # convergence
        if it >= min_iter and _logLik_new - _logLik_old < epsilon_conv:
            if verbose:
                info("EM finished in %d iterations with logLik %.4e" 
                      %(it, _logLik_new))
            break
        elif _logLik_new < _logLik_old:
            if verbose:
                warn("likelihood decreases in EM algorithm!")
                info("    %d %f %f" % (it, _logLik_old, _logLik_new))
    
    # soft phasing bins counts: `ad_sum`
    ad_sum = Z.T @ AD + (1 - Z.T) @ BD

    # hard phasing bins counts:`ad_sum1``
    Z_argmax = np.argmax(Z, axis = 1)
    Z_hard_assign = np.vstack((1 - Z_argmax, Z_argmax)).T
    ad_sum1 = Z_hard_assign.T @ AD + (1 - Z_hard_assign.T) @ BD 
    dp_sum = DP.sum(axis = 0)
    
    return ad_sum, ad_sum1, dp_sum, Z, thetas, _logLik_new



def gaussian_smoothing_1d(v, x, b, a = 0):
    """1-D smoothing with gaussian kernel.

    The kernel is exp(a - (x1-x2)^2 / b^2).

    Parameters
    ----------
    v : array-like
        The values to be smoothed.
    x : array-like
        The 1-d coordinates, same length with `v`.
    b : float
        The kernel parameter controling the spread of the Gaussian.
    a : float
        An scaling coefficient.
        The kernel can be reformulated as exp(a)exp((x1-x2)^2 / b^2).

    Returns
    -------
    array-like
        The smoothed values.
    """
    u = v.copy()
    for i in range(len(v)):
        w = np.exp(a - (x - x[i])**2 / b**2)
        w = w / np.sum(w)
        u[i] = np.sum(v * w)
    return(u)



def loglik_amplify(X, axis = -1):
    """
    Amplify the log likelihood matrix by subtract the maximum.
    X should be numpy.array, otherwise will be transformed to it.
    
    Example
    -------
    X = np.random.rand(3, 5, 8)
    loglik_amplify(X, axis = 1)
    """
    shape2 = list(X.shape)
    shape2[axis] = 1
    
    if type(X) == np.matrix:
        X = X.toarray()
    
    X_max = np.max(X, axis = axis, keepdims = True)
    return X - X_max



def normalize(X, axis = -1):
    """
    Normalization of tensor with sum to 1.
    X should be numpy.array, otherwise will be transformed to it.
    
    Example
    -------
    X = np.random.rand(3, 5, 8)
    tensor_normalize(X, axis = 1)
    """
    shape2 = list(X.shape)
    shape2[axis] = 1
    
    if type(X) == np.matrix:
        X = X.toarray()
    if sp.sparse.issparse(X):
        X = X.toarray()
    
    X_sum = np.sum(X, axis = axis, keepdims = True)
    return X / X_sum



# the default value of `b` is estimated in a simulated ST dataset.
def snp_smooth_gk(z, c, b = 20000, a = 0):
    """SNP haplotype assignment smoothing with gaussian kernel.
    
    Parameters
    ----------
    z : array-like
        The probabilities of haplotype assignment to be smoothed.
    c : array-like
        The SNP genomic coordinates, same length with `z`.
    b : float
        The kernel parameter controling the spread of the Gaussian.
    a : float
        An scaling coefficient.
        The kernel can be reformulated as exp(a)exp((x1-x2)^2 / b^2).

    Returns
    -------
    array-like
        The smoothed values.
    """
    return gaussian_smoothing_1d(z, c, b = b, a = a)
