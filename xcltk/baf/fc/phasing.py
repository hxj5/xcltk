# phasing.py - local phasing in regions.
# further check the SNP-level reference phasing states after running Eagle2.


import numpy as np
import scipy as sp
from logging import info
from logging import warning as warn
from ..localphase import snp_local_phasing


    
def reg_local_phasing(
    reg,
    AD, DP,
    kws_localphase = None,
    verbose = False
):
    if kws_localphase is None:
        kws_localphase = dict()

    cell_idx = DP.sum(axis = 1) > 0
    snp_idx = DP.sum(axis = 0) > 0
    AD, DP = AD[np.ix_(cell_idx, snp_idx)], DP[np.ix_(cell_idx, snp_idx)]
    reg.snp_list = [s for i, s in zip(snp_idx, reg.snp_list) if i]

    BD = DP - AD
    flip = np.array([snp.ref_idx == 1 for snp in reg.snp_list])
    AD_ref_phased = AD * (1 - flip.T) + BD * flip.T

    res = snp_local_phasing(
        AD = AD_ref_phased, 
        DP = DP,
        positions = np.array([s.pos for s in reg.snp_list]),
        verbose = verbose, 
        **kws_localphase
    )
    if res is None:
        warn("local phasing for region '%s' failed!" % reg.name)
        return(reg)
    
    flip, ad_sum, ad_sum1, dp_sum, Z, thetas, logLik_new = res
    if np.mean(flip) > 0.5:       # try to flip smaller number of SNPs.
        flip = 1 - flip
    else:
        flip = flip + 0
    flip = flip.astype(int)
    
    assert len(reg.snp_list) == len(flip)
    for i in range(len(reg.snp_list)):
        if flip[i] == 1:
            snp = reg.snp_list[i]
            snp.ref_idx, snp.alt_idx = 1 - snp.ref_idx, 1 - snp.alt_idx
            snp.gt = {snp.ref:snp.ref_idx, snp.alt:snp.alt_idx}
    return(reg)
