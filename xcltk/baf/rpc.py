# rpc.py - reference phasing correction


import anndata as ad
import getopt
import multiprocessing
import numpy as np
import os
import pandas as pd
import scipy as sp
import scipy.io
import sys
import time

from scipy.special import logsumexp
from ..config import APP, VERSION



COMMAND = "rpc"


def format_theta(theta, epsilon = 1e-6, inplace = True):
    """Format Theta Values
    
    Add a tiny value to theta = 0 or subtract a tiny value from theta = 1, to
    avoid numeric issues in calculating log(theta) or log(1 - theta).

    Parameters
    ----------
    theta: A float matrix containing 2 columns. Theta values, the sum
        of each row is 1 [np.ndarray((M, 2))]
    epsilon: A float scalar. The tiny value added or subtracted [float]
    inplace: A bool. Whether to modify `theta` in-place [bool]

    Returns
    -------
    A formatted float matrix of theta(s), in the same dimension with
    the input matrix [np.ndarray((M, 2))]

    Example
    -------
    theta = np.zeros((10, 2))
    theta[:, 0] = np.random.rand(10)
    theta[0, 0] = 0
    theta[:, 1] = 1 - theta[:, 0]

    theta_new = format_theta(theta, epsilon = 1e-10)
    print(theta_new)
    """
    if not inplace:
        theta = theta.copy()

    idx = (theta[:, 0] <= 0) + (theta[:, 1] >= 1)
    theta[idx, 0] = epsilon
    theta[idx, 1] = 1 - epsilon
    
    idx = (theta[:, 0] >= 1) + (theta[:, 1] <= 0)
    theta[idx, 0] = 1 - epsilon
    theta[idx, 1] = epsilon
    
    return(theta)



def ref_phasing_correction(A, D, hap, epsilon_hap = 1e-2, 
                           max_iter = 1000, epsilon_conv = 1e-6,
                           verbose = True):
    """XClone Reference Phasing Correction
    
    Correct errors in XClone reference phasing, taking AD, DP matrix, and
    haplotypes as input.
    
    Parameters
    ----------
    A: An integer AD count matrix of N blocks and M cells [np.ndarray((N, M))]
    D: An integer DP count matrix of N blocks and M cells [np.ndarray((N, M))]
    hap: A binary vector. Haplotype (0 for H0, 1 for H1) of ALT alleles in 
        N blocks before correction [np.array((N, ))]
    epsilon_hap: A float. Estimated error rate of reference phased haplotypes,
        used in constructing prior `pi`. Note that `pi` is a float matrix for
        N blocks [np.ndarry((N, 2))], which are the prior probabilities of 
        the ALT allele being assigned to H0 or H1 in each block [float]
    max_iter: An integer scalar. Maximum iterations of EM [int]
    epsilon_conv: A float scalar. The minimum absolute difference between
        two EM iterations required to stop iterations [float]
    verbose: A bool. Whether to output detailed logging information [bool]

    Returns
    -------
    A dict of 7 elements. 
    psi: an N x 2 float matrix, the expectations of the latent variables, 
        each row is the estimated posterial probabilities of the ALT allele 
        being assigned to H0 or H1 [np.ndarry((N, 2))]
    theta: a M x 2 float matrix, each row is the probability parameters 
        in binomial distributions for H0 and H1 [np.ndarry((M, 2))]
    pi: a float matrix, which are the prior probabilities of the ALT allele 
        being assigned to H0 or H1 in each block [np.ndarry((N, 2))]
    log_lik_list: an float vector containing loglikelihood from each 
        iteration [np.ndarray((n_iter, ))]
    flip: an integer vector of N elements indicating whether AD in each 
        block should be flipped [np.ndarray((N, 1))]
    AD_phased: an integer phased/flipped AD count matrix of N blocks and 
        M cells [np.ndarray((N, M))]
    hap_new: A binary vector. Haplotype (0 for H0, 1 for H1) of ALT alleles in
        N blocks after correction [np.array((N, ))]
        
    Example
    -------
    dat = run_simulation(N = 100, M = 220, pi = 0.3, DP_lambda = 1226,
                         perc = 0.8, 
                         theta_mu2 = c(0.45, 0.9), theta_sigma = 0.15)
    res = ref_phasing_correction(A = dat["AD"], D = dat["DP"], 
                                 hap = dat["hap"])
    print(res["hap_new"])
    """
    func = "ref_phasing_correction"
    tiny = 1e-10
  
    N, M = A.shape
    B = D - A
    
    colsum_D = np.array(D.sum(axis = 0).reshape((-1, 1)), dtype = "float")
    colsum_D[colsum_D <= 0] = tiny
  
    # Initialization
    theta = np.zeros((M, 2))
    theta[:, 0] = np.random.rand(M)
    theta[:, 1] = 1 - theta[:, 0]
    theta = format_theta(theta, tiny)
  
    pi = np.zeros((N, 2))
    pi[:, 0] = [1 - epsilon_hap if h == 0 else epsilon_hap for h in hap]
    pi[:, 1] = 1 - pi[:, 0]
  
    phi = A @ np.log(theta) + B @ np.log(1 - theta) + np.log(pi)
    xi = logsumexp(phi, axis = 1)
    log_lik = np.sum(xi)
  
    log_lik_list = [log_lik]
    if verbose:
        print("[%s][init] log-likelihood = %f." % (func, log_lik))

    itr = 0
    while itr < max_iter:
        # E-step
        psi = np.exp(phi - xi.reshape((-1, 1)))
    
        # M-step
        theta = (A.T @ psi + B.T @ (1 - psi)) / colsum_D
        theta = format_theta(theta, tiny)
        
        # Check
        phi = A @ np.log(theta) + B @ np.log(1 - theta) + np.log(pi)
        xi = logsumexp(phi, axis = 1)
        log_lik_new = np.sum(xi)
    
        log_lik_list.append(log_lik_new)
        if verbose:
            print("[%s][iter-%d] log-likelihood = %f." %
                  (func, itr, log_lik_new))
    
        if abs(log_lik - log_lik_new) < epsilon_conv:
            if verbose:
                print("[%s][stop] abs(log_lik - log_lik_new) < epsilon (%e)." %
                      (func, epsilon_conv))
                print("[%s][quit] final log-likelihood = %f, iter = %d." %
                       (func, log_lik_new, itr))
            break
        else:
            log_lik = log_lik_new
        itr = itr + 1
  
    if verbose and itr >= max_iter:
        print("[%s][stop] number of iterations exceeds max_iter (%d)." % 
               (func, max_iter))
        print("[%s][quit] final log-likelihood = %f, iter = %d." %
               (func, log_lik_new, itr))
  
    hap_new = np.argmax(psi, axis = 1)
    flip = 1 - (hap == hap_new).reshape((-1, 1))
    AD_phased = A
  
    return(dict(
        psi = psi, theta = theta, pi = pi,
        log_lik_list = log_lik_list,
        flip = flip, AD_phased = AD_phased, hap_new = hap_new))



def multi_init_rpc(A, D, hap, epsilon_hap = 1e-2, 
                   max_iter = 1000, epsilon_conv = 1e-6,
                   n_init = 100, verbose = True):
    """XClone Reference Phasing Correction with Multiple Random Initialization
    
    Correct errors in XClone reference phasing with multiple random 
    initialization, taking AD, DP matrices, and hap as input, to alleviate 
    potential issue of local optima.

    Parameters
    ----------
    @inheritParams ref_phasing_correction
    n_init: An integer scalar. Number of random initialization [int]
    
    Returns
    -------
    A dict of 3 elements.
    idx_max: An integer. The index (0-based) of the maximum loglikelihood
        in the list `all_log_lik` [int]
    all_log_lik: a numeric vector of final loglikelihood for each random 
        initialization [np.ndarray((n_init, ))]
    all_res: a list of n_init elements, each element is the returned value
        of `run_local_phasing` for corresponding initialization [list]

    Example
    -------
    dat = run_simulation(N = 100, M = 220, pi = 0.3, DP_lambda = 1226,
                         perc = 0.8, 
                         theta_mu2 = c(0.45, 0.9), theta_sigma = 0.15)
    res = multi_init_rpc(A = dat["AD"], D = dat["DP"], hap = dat["hap"], 
                         n_init = 100)
    """
    func = "multi_init_rpc"
  
    all_res = []
    all_log_lik = np.zeros((n_init, ))
    for idx in range(n_init):
        res = ref_phasing_correction(
            A = A, D = D, hap = hap, epsilon_hap = epsilon_hap,
            max_iter = max_iter, epsilon_conv = epsilon_conv, verbose = False)
        final_log_lik = res["log_lik_list"][-1]
        if verbose:
            print("[%s][init-%d] final loglikelihood is %f." %
                   (func, idx, final_log_lik))
        all_res.append(res)
        all_log_lik[idx] = final_log_lik
    idx_max = np.argmax(all_log_lik)
    if verbose:
        print("[%s][stop] max loglikelihood is %f from initialization-%d." %
               (func, np.max(all_log_lik), idx_max))
    return(dict(idx_max = idx_max, all_log_lik = all_log_lik, 
                all_res = all_res))



def region_rpc(adata,
               reg_name, reg_chrom, reg_start, reg_end,
               n_init = 100, epsilon_hap = 1e-2,
               verbose = False):
    """Reference Phasing Correctioin in Specific Region

    Parameters
    ----------
    adata: AnnData object storing SNP x Cell data. Its `obs` contains 6 columns
        chrom, pos, ref, alt, ref_hap, alt_hap obtained from calling 
        `load_snp_meta` [ad.AnnData]
    reg_name: region name/ID [str]
    reg_chrom: region chromosome [str]
    reg_start: 1-based start position of region, inclusive [int]
    reg_end: 1-based end position of region, inclusive [int]
    n_init: number of random initialization for EM algorithm [int]
    epsilon_hap: A float. Estimated error rate of reference phased haplotypes,
        used in constructing prior `pi` [float]
    verbose: A bool. Whether to output detailed logging information [bool]

    Returns
    -------
    A dict containing 7 elements [dict]:
    reg_name, reg_chrom, reg_start, reg_end, as demonstrated in `Parameters`.
    n_flip: number of SNPs whose haplotypes should be flipped [int]
    n_snp: number of total SNPs [int]
    snp_meta: meta information of all SNPs. A dataframe containing 8 columns:
        chrom, chromosome name [str]
        pos, 1-based genomic position of SNP [int]
        ref, REF allele, one of A/C/G/T/N [str]
        alt, ALT allele, one of A/C/G/T/N [str]
        ref_hap, haplotype ID of REF allele [int]
        alt_hap, haplotype ID of ALT allele [int]
        ref_hap_new, corrected haplotype ID of REF allele [int]
        alt_hap_new, corrected haplotype ID of ALT allele [int]
    """
    func = "region_rpc"
    if verbose:
        sys.stderr.write("[I::%s] processing region %s ...\n" % \
            (func, reg_name))

    reg_dat = adata[(adata.obs.chrom == reg_chrom) &
                    (adata.obs.pos >= reg_start) &
                    (adata.obs.pos <= reg_end)]

    AD = reg_dat.layers["AD"].toarray()
    DP = reg_dat.layers["DP"].toarray()
    hap = np.array(reg_dat.obs.alt_hap, dtype = "int")
    N, M = AD.shape

    if N <= 0:
        n_flip = 0
        snp_meta = pd.DataFrame([], columns = ["chrom", "pos", "ref", "alt", \
            "ref_hap", "alt_hap", "ref_hap_new", "alt_hap_new"])
    elif N == 1:
        n_flip = 0
        snp_meta = reg_dat.obs.copy()
        snp_meta["ref_hap_new"] = 1 - hap
        snp_meta["alt_hap_new"] = hap
    else:
        # correction
        mi_res = multi_init_rpc(
            AD, DP, hap, epsilon_hap = epsilon_hap, 
            max_iter = 1000, epsilon_conv = 1e-6,
            n_init = n_init, verbose = False)

        res = mi_res["all_res"][mi_res["idx_max"]]
        snp_meta = reg_dat.obs.copy()
        snp_meta["ref_hap_new"] = 1 - res["hap_new"]
        snp_meta["alt_hap_new"] = res["hap_new"]    
        n_flip = np.sum(res["flip"])

    if verbose:
        sys.stderr.write(
            "[I::%s] region %s: %d (total %d) SNPs are flipped.\n" % \
            (func, reg_name, n_flip, N))

    return(dict(reg_name = reg_name, reg_chrom = reg_chrom,
                reg_start = reg_start, reg_end = reg_end,
                n_flip = n_flip, n_snp = N, 
                snp_meta = snp_meta))



def show_progress(res = None):
    return(res)



# TODO: use region.RegionSet to implement @param gene_set.
def multi_reg_rpc(
    adata, reg_set,
    n_init = 100, epsilon_hap = 1e-2,
    n_proc = 1, verbose = False):
    """Reference Phasing Correctioin in Multiple Regions

    Parameters
    ----------
    adata: AnnData object storing SNP x Cell data. Its `obs` contains 6 columns
        chrom, pos, ref, alt, ref_hap, alt_hap obtained from calling 
        `load_snp_meta` [ad.AnnData]
    reg_set: a dataframe containing region meta information. It has 4 columns
        chrom, start, end, name obtained from calling `load_region_meta`
        [pd.DataFrame]
    n_init: number of random initialization for EM algorithm [int]
    epsilon_hap: A float. Estimated error rate of reference phased haplotypes,
        used in constructing prior `pi` [float]
    n_proc: number of processes [int]
    verbose: A bool. Whether to output detailed logging information [bool]

    Returns
    -------
    A dict of 2 elements:
    snp_stat, a dataframe containing 8 columns [pd.DataFrame]:
        chrom, pos, ref, alt, ref_hap, alt_hap, ref_hap_new, alt_hap_new,
        see `region_rpc` for details.
    reg_stat, a dataframe containing 6 columns [pd.DataFrame]:
        chrom [str], start [int], end [int], name [str] of region 
        n_flip: number of SNPs whose haplotypes should be flipped [int]
        n_snp: number of total SNPs [int]
    """
    func = "multi_reg_rpc"

    n_reg = reg_set.shape[0]
    result = []
    if n_proc > 1:
        pool = multiprocessing.Pool(processes = n_proc)
        for i in range(n_reg):
            res = pool.apply_async(region_rpc, (
                adata,                          
                reg_set["name"][i], reg_set["chrom"][i],        
                reg_set["start"][i], reg_set["end"][i],
                n_init,
                epsilon_hap,
                verbose), 
                callback = show_progress)
            result.append(res)
        pool.close()
        pool.join()
        result = [res.get() for res in result]
    else:
        for i in range(n_reg):
            res = region_rpc(
                adata,
                reg_set["name"][i], reg_set["chrom"][i],
                reg_set["start"][i], reg_set["end"][i],
                n_init,
                epsilon_hap,
                verbose = verbose)
            result.append(res)

    assert n_reg == len(result)

    reg_stat, snp_stat = None, None

    reg_stat = pd.DataFrame([res["reg_chrom"] for res in result], 
                            columns = ["chrom"])
    reg_stat["start"] = np.array([res["reg_start"] for res in result])
    reg_stat["end"] = np.array([res["reg_end"] for res in result])
    reg_stat["name"] = np.array([res["reg_name"] for res in result])
    reg_stat["n_flip"] = np.array([res["n_flip"] for res in result])
    reg_stat["n_snp"] = np.array([res["n_snp"] for res in result])

    reg_stat = reg_stat.sort_values(by = ["chrom", "start", "end", "name"])

    snp_stat = pd.DataFrame([], columns = ["chrom", "pos", "ref", "alt", \
        "ref_hap", "alt_hap", "ref_hap_new", "alt_hap_new"])
    for i in range(n_reg):
        snp_stat = pd.concat([snp_stat, result[i]["snp_meta"]], axis = 0)

    snp_stat = snp_stat.sort_values(by = ["chrom", "pos"])

    return((snp_stat, reg_stat))



class Config:
    def __init__(self):
        self.def_epsilon_hap = 1e-2
        self.def_n_init = 100
        self.def_n_proc = 1

        self.barcode_fn = None
        self.AD_mtx_fn = None
        self.DP_mtx_fn = None
        self.out_dir = None
        self.snp_fn = None
        self.region_fn = None
        self.sid = None
        self.verbose = False

        self.epsilon_hap = self.def_epsilon_hap
        self.n_init = self.def_n_init
        self.n_proc = self.def_n_proc


        
def load_cell_meta(barcode_fn):
    """Load Cell Meta Information

    Parameters
    ----------
    barcode_fn: path to cell-barcode file without header [str]

    Returns
    -------
    A dataframe containing 1 column [pd.DataFrame]: 
    cell, cell barcodes [str]
    """
    cell_meta = pd.read_csv(barcode_fn, sep = "\t", header = None)
    cell_meta.columns = ["cell"]
    return(cell_meta)



def load_region_meta(region_fn):
    """Load Region Meta Information

    Parameters
    ----------
    region_fn: path to TSV file that containing 4 columns without header:
        chrom, chromosome name [str]
        start, 1-based start position of region, inclusive [int]
        end, 1-based end position of region, inclusive [int]
        name, region name or ID [str]

    Returns
    -------
    A dataframe containing 4 column [pd.DataFrame]:
        chrom [str], start [int], end [int], name [str]
    """
    region_meta = pd.read_csv(region_fn, sep = "\t", header = None, 
                            dtype = {0: str})
    region_meta.columns = ["chrom", "start", "end", "name"]
    region_meta["chrom"] = region_meta["chrom"].astype(str)
    region_meta["chrom"] = region_meta["chrom"].str.replace(  \
        "^chr", "", regex = True)
    return(region_meta)



def load_snp_meta(snp_fn):
    """Load SNP Meta Information

    Parameters
    ----------
    snp_fn: path to VCF file that containing the phased GT [str]

    Returns
    -------
    A dataframe containing 6 columns [pd.DataFrame]: 
    chrom, chromosome name [str]
    pos, 1-based genomic position of SNP [int]
    ref, REF allele, one of A/C/G/T/N [str]
    alt, ALT allele, one of A/C/G/T/N [str]
    ref_hap, haplotype ID of REF allele [int]
    alt_hap, haplotype ID of ALT allele [int]
    """
    snp_meta = pd.read_csv(snp_fn, sep = "\t", comment = "#", header = None,
                            dtype = {0: str})
    snp_meta = snp_meta[[0, 1, 3, 4, 9]]
    snp_meta.columns = ["chrom", "pos", "ref", "alt", "gt"]
    snp_meta["chrom"] = snp_meta["chrom"].astype(str)
    snp_meta["chrom"] = snp_meta["chrom"].str.replace("^chr", "", regex = True)
    snp_meta["gt"] = snp_meta["gt"].str.split(":").str[0]
    snp_meta["ref_hap"] = snp_meta["gt"].str.split("|").str[0]
    snp_meta["ref_hap"] = np.array(snp_meta["ref_hap"], dtype = "int")
    snp_meta["alt_hap"] = snp_meta["gt"].str.split("|").str[1]
    snp_meta["alt_hap"] = np.array(snp_meta["alt_hap"], dtype = "int")
    snp_meta = snp_meta.drop("gt", axis = 1)
    return(snp_meta)



def load_AD(AD_mtx_fn):
    AD = sp.io.mmread(AD_mtx_fn).tocsr()
    return(AD)


def load_DP(DP_mtx_fn):
    DP = sp.io.mmread(DP_mtx_fn).tocsr()
    return(DP)



def usage(conf, fp = sys.stdout):
    s =  "\n"
    s += "Version: %s\n" % VERSION
    s += "Usage:   %s %s <options>\n" % (APP, COMMAND)
    s += "\n" 
    s += "Options:\n"
    s += "  -b, --barcode FILE     A plain file listing all effective cell barcode.\n"
    s += "  -p, --nproc INT        Number of processes [%d]\n" % conf.def_n_proc
    s += "  -A, --AD FILE          AD sparse matrix file.\n"
    s += "  -D, --DP FILE          DP sparse matrix file.\n"
    s += "  -O, --outdir DIR       Output directory.\n"
    s == "  -P, --phasedSNP FILE   A VCF file listing phased SNPs (i.e., containing phased GT).\n"
    s += "  -R, --region FILE      A TSV file with 4 columns <chrom> <start> <end> <name>.\n"
    s += "  -S, --sid STR          Sample ID.\n"
    s += "  -h, --help             Print this message and exit.\n"
    s += "  -v, --verbose          If set, print detailed logging information.\n"
    s += "  --epsilonHAP FLOAT     Estimated error rate of reference phased haplotypes [%e]\n" % conf.def_epsilon_hap
    s += "  --ninit INT            Number of random initialization [%d]\n" % conf.def_n_init
    s += "\n"

    fp.write(s)

    

def check_args(conf):
    func = "check_args"

    if conf.barcode_fn is None:
        sys.stderr.write("[E::%s] barcode file is needed!\n" % func)
        return(-1)

    if conf.AD_mtx_fn is None:
        sys.stderr.write("[E::%s] AD matrix file is needed!\n" % func)
        return(-1)

    if conf.DP_mtx_fn is None:
        sys.stderr.write("[E::%s] DP matrix file is needed!\n" % func)
        return(-1)

    if conf.out_dir is None:
        sys.stderr.write("[E::%s] output dir is needed!\n" % func)
        return(-1)
    elif not os.path.isdir(conf.out_dir):
        os.mkdir(conf.out_dir)

    if conf.snp_fn is None:
        sys.stderr.write("[E::%s] SNP file is needed!\n" % func)
        return(-1)

    if conf.region_fn is None:
        sys.stderr.write("[E::%s] region file is needed!\n" % func)
        return(-1)

    return(0)



def quit(main_func, cmdline, start_time, state = 0, err_msg = None):
    func = main_func

    if state >= 0:
        sys.stdout.write("[I::%s] All Done!\n" % func)
    else:
        sys.stderr.write("[E::%s] %s\n" % (func, err_msg))
        sys.stdout.write("[E::%s] Running program failed.\n" % func)
        sys.stdout.write("[E::%s] Quiting ...\n" % func)

    sys.stdout.write("[I::%s] CMD: %s\n" % (func, cmdline))

    end_time = time.time()
    time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(end_time))
    sys.stdout.write("[I::%s] end time: %s\n" % (func, time_str))
    sys.stdout.write("[I::%s] time spent: %.2fs\n" % (func, end_time - start_time))
    sys.stdout.write("\n")


    
def rpc_main(argv):
    func = "rpc_main"

    start_time = time.time()

    conf = Config()
    if len(argv) <= 2:
        usage(conf, sys.stdout)
        sys.exit(0)

    opts, args = getopt.getopt(argv[2:], 
        "-b:-p:-A:-D:-O:-P:-R:-S:-h-v", 
        [
            "barcode=", "nproc=", 
            "AD=", "DP=",
            "outdir=", "phasedSNP=", 
            "region=",  "sid=", "help", "verbose",
            "epsilonHAP=", "ninit="
        ])

    for op, val in opts:
        if len(op) > 2:
            op = op.lower()
        if op in   ("-b", "--barcode"): conf.barcode_fn = val
        elif op in ("-p", "--nproc"): conf.n_proc = int(val)
        elif op in ("-A", "--ad"): conf.AD_mtx_fn = val
        elif op in ("-D", "--dp"): conf.DP_mtx_fn = val
        elif op in ("-O", "--outdir"): conf.out_dir = val
        elif op in ("-P", "--phasedsnp"): conf.snp_fn = val
        elif op in ("-R", "--region"): conf.region_fn = val
        elif op in ("-S", "--sid"): conf.sid = val
        elif op in ("-h", "--help"): usage(conf); sys.exit(0)
        elif op in ("-v", "--verbose"): conf.verbose = True
        elif op in ("--epsilonhap"): conf.epsilon_hap = float(val)
        elif op in ("--ninit"): conf.n_init = int(val)

    sys.stdout.write("\n")
    time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(start_time))
    sys.stdout.write("[I::%s] start time: %s.\n" % (func, time_str))

    cmdline = " ".join(argv)
    sys.stdout.write("[I::%s] CMD: %s\n" % (func, cmdline))

    # check args
    sys.stdout.write("[I::%s] check args ...\n" % func)

    if check_args(conf) < 0:
        quit(func, cmdline, start_time, state = -3, 
             err_msg = "check_args failed!")
        return(-3)

    # load data
    sys.stdout.write("[I::%s] load data ...\n" % func)

    AD = load_AD(conf.AD_mtx_fn)
    DP = load_DP(conf.DP_mtx_fn)
    cell_meta = load_cell_meta(conf.barcode_fn)
    snp_meta = load_snp_meta(conf.snp_fn)
    reg_set = load_region_meta(conf.region_fn)

    adata = ad.AnnData(AD, obs = snp_meta, var = cell_meta)
    adata.layers["AD"] = AD
    adata.layers["DP"] = DP

    n_snp, n_cell = AD.shape
    n_reg = reg_set.shape[0]
    sys.stdout.write("[I::%s] %d SNPs, %d cells, %d regions loaded.\n" % \
                     (func, n_snp, n_cell, n_reg))

    # reference phasing correction
    sys.stdout.write("[I::%s] reference phasing correction ...\n" % func)

    snp_stat, reg_stat = multi_reg_rpc(
            adata, reg_set,
            n_init = conf.n_init, epsilon_hap = conf.epsilon_hap,
            n_proc = conf.n_proc, verbose = conf.verbose         
        )

    if snp_stat is None or reg_stat is None:
        quit(func, cmdline, start_time, state = -5,
             err_msg = "reference phasing correction for multiple regions failed!")
        return(-5)

    # output
    sys.stdout.write("[I::%s] output results ...\n" % func)

    out_fn_snp = os.path.join(conf.out_dir, "%s_snp_stat.tsv" % conf.sid)
    out_fn_reg = os.path.join(conf.out_dir, "%s_region_stat.tsv" % conf.sid)
    snp_stat.to_csv(out_fn_snp, sep = "\t", index = False)
    reg_stat.to_csv(out_fn_reg, sep = "\t", index = False)
    
    sys.stdout.write("[I::%s] results saved to dir '%s'.\n" % \
        (func, conf.out_dir))

    # quit
    quit(func, cmdline, start_time, state = 0)
    return(0)
