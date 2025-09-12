# refphase.py - reference phasing.


import multiprocessing
import os
import stat
import subprocess
import sys

from logging import error, info
from ..utils.base import assert_e, assert_n



def ref_phasing(
    target_vcf_list, ref_vcf_list, out_prefix_list,
    gmap_fn,
    eagle_fn,
    out_dir,
    ncores = 1,
    script_fn_prefix = None, log_fn_prefix = None,
    verbose = False
):
    """Reference phasing

    The functions takes as input the target and reference VCF files, together
    with the reference panels for SNP phasing, and outputs corresponding
    phased VCF (.vcf.gz) files in BGZF format.
    It internally calls the `Eagle2` command-line tool for phasing.

    Parameters
    ----------
    target_vcf_list : list
        A list of target VCF files, i.e., the input VCF files to be phased.
    ref_vcf_list : list
        A list of reference VCF files, i.e., the input VCF files containing
        reference genotypes.
    out_prefix_list : list
        A list of prefixes of the output VCF files. The output phased VCF
        file path will be, e.g., `<out_prefix_list[0]>.vcf.gz`.
    gmap_fn : str
        The genetic map file.
    eagle_fn : str
        The Eagle2 binary executable file.
    out_dir : str
        The output dir.
    ncores : int
        Number of threads.
    script_fn_prefix : str
        Prefix to the script files that run Eagle2. If `None`, use default
        value `<out_dir>/run_phasing`.
    log_fn_prefix : str
        Prefix to the logging files that record the output of Eagle2. 
        If `None`, use default value `<out_dir>/phasing`.
    verbose : bool
        Whether to show detailed logging information.

    Returns
    -------
    None

    Raises
    ------
    ValueError
        When some input arguments are invalid.
    AssertionError
        When some input arguments are invalid.
    """
    if verbose:
        info("start reference phasing ...")

    # check args
    if not target_vcf_list:
        raise ValueError("target vcf list is empty.")
    if not ref_vcf_list:
        raise ValueError("reference vcf list is empty.")
    if not out_prefix_list:
        raise ValueError("output prefix list is empty.")
    if len(target_vcf_list) != len(ref_vcf_list):
        raise ValueError("number of target and reference vcf files should be the same.")
    if len(target_vcf_list) != len(out_prefix_list):
        raise ValueError("number of target vcf files and output prefix should be the same.")
    for fn in target_vcf_list:
        assert_e(fn)
    for fn in ref_vcf_list:
        assert_e(fn)
    for prefix in out_prefix_list:
        assert_n(prefix)

    assert_e(gmap_fn)
    assert_e(eagle_fn)
    assert_n(out_dir)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok = True)

    if script_fn_prefix is None:
        script_fn_prefix = os.path.join(out_dir, "run_phasing")
    if log_fn_prefix is None:
        log_fn_prefix = os.path.join(out_dir, "phasing")

        
    # run phasing
    if verbose:
        info("run phasing ...")

    pool = multiprocessing.Pool(processes = ncores)
    result = []
    idx = 0
    for target_vcf_fn, ref_vcf_fn, out_prefix in zip(
        target_vcf_list, ref_vcf_list, out_prefix_list
    ):
        cmd  = ""
        cmd += "%s  \\\n" % eagle_fn
        cmd += "    --vcfTarget  %s    \\\n" % target_vcf_fn 
        cmd += "    --vcfRef  %s   \\\n" % ref_vcf_fn 
        cmd += "    --geneticMapFile  %s    \\\n" % gmap_fn
        cmd += "    --outPrefix  %s    \\\n" % out_prefix
        cmd += "    --numThreads  %d    \n" % 1        # min(ncores, 3)
        cmd += "\n"

        idx += 1
        script_fn = script_fn_prefix + "_%d.sh" % idx
        log_fn = log_fn_prefix + "_%d.log" % idx
        with open(script_fn, "w") as fp:
            fp.write(cmd)
        st = os.stat(script_fn)
        os.chmod(script_fn, st.st_mode | stat.S_IXUSR)

        result.append(pool.apply_async(
            func = ref_phasing1, 
            args = (script_fn, log_fn))
        )
    pool.close()
    pool.join()

    

def ref_phasing1(script_fn, log_fn):
    ret = None
    try:
        proc = subprocess.Popen(
            args = "%s 2>&1 | tee %s" % (script_fn, log_fn),
            shell = True,
            executable = "/bin/bash", 
            stdout = subprocess.PIPE, 
            stderr = subprocess.PIPE
        )
        outs, errs = proc.communicate()
        ret = proc.returncode
        if ret != 0:
            raise RuntimeError(str(errs.decode()))
    except Exception as e:
        error(str(e))
        error("Error: phasing failed (retcode '%s')." % str(ret))
        sys.exit(1)
