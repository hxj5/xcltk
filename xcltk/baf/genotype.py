# genotype.py - preprocess the input BAM file to generate reference-phased cell x gene AD & DP matrices.

import gzip
import os
import pandas as pd
import random
import stat
import subprocess
import sys

from logging import error, info
from logging import warning as warn
from ..utils.base import assert_e, assert_n


def pileup(
    sam_fn = None, sam_list_fn = None, 
    barcode_fn = None, sample = None,
    snp_vcf_fn = None,
    out_dir = None,
    mode = "10x",
    cell_tag = "CB", umi_tag = "UB",
    ncores = 1,
    min_maf = 0.1, min_count = 20,
    script_fn = None, log_fn = None
):
    """Pileup indexed BAM file, supporting both single-cell and bulk data.

    The function internally will call `cellsnp-lite` for pileup.
    The input BAM file(s) should be specified by one and only one of `sam_fn`
    and `sam_list_fn`.

    Parameters
    ----------
    sam_fn : str
        Comma separated indexed BAM/CRAM file(s).
    sam_list_fn : str
        A file listing BAM/CRAM files, each per line.
    barcode_fn : str
        A plain file listing all effective cell barcodes (for 10x) or cell IDs (for smartseq).
    sample : str
        Sample ID (for bulk data).
    snp_vcf_fn : str
        A vcf file listing all candidate SNPs.
    out_dir : str
        Output dir.
    mode : str
        One of "10x", "smartseq", and `bulk`.
    cell_tag : str
        Cell barcode tag, set to `None` to turn it off.
    umi_tag : str
        UMI tag, set to `None` to turn it off.
    ncores : int
        Number of threads.
    min_maf : float
        Minimum minor allele frequency.
    min_count : int
        Minimum aggregated count.
    script_fn : str
        Path to the script file that runs cellsnp-lite. If `None`, use default
        path `<out_dir>/run_pileup.sh`.
    log_fn : str
        Path to the logging file that records the output of cellsnp-lite. If
        `None`, use default path `<out_dir>/pileup.log`.

    Returns
    -------
    None

    Raises
    ------
    AssertionError
        when some input file or dirs are invalid.

    Notes
    -----
    1. For smartseq data, the order of the BAM files (in `sam_fn` or `sam_list_fn`)
       and the sample IDs (in `barcode_fn`) should match each other.
    """
    # check args
    assert mode in ("10x", "smartseq", "bulk")

    sam_list = None
    if sam_fn is None:
        assert_e(sam_list_fn)
        sam_list = [line.strip() for line in open(sam_list_fn, "r")]
        for fn in sam_list:
            assert_e(fn)
    else:
        assert sam_list_fn is None
        sam_list = sam_fn.split(",")
        for fn in sam_list:
            assert_e(fn)
    assert len(sam_list) >= 1

    if mode == "bulk":
        assert_n(sample)
        assert len(sam_list) == 1
    else:
        assert_e(barcode_fn)
        barcodes = [line.strip() for line in open(barcode_fn, "r")]
        if mode == "smartseq":
            assert len(sam_list) == len(barcodes)

    assert_e(snp_vcf_fn)

    if mode in ("smartseq", "bulk"):
        if str(cell_tag) != "None":
            warn("cell tag is not 'None' in mode '%s'." % mode)
        if str(umi_tag) != "None":
            warn("umi tag is not 'None' in mode '%s'." % mode)
    else:
        if str(cell_tag) == "None":
            warn("cell tag is 'None' in mode '%s'." % mode)
        if str(umi_tag) == "None":
            warn("umi tag is 'None' in mode '%s'." % mode)

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    
    if not script_fn:
        script_fn = os.path.join(out_dir, "run_pileup.sh")
    if not log_fn:
        log_fn = os.path.join(out_dir, "pileup.log")

    # generate pileup script
    cmd  = "cellsnp-lite  \\\n"
    if sam_fn:
        cmd += "    -s  %s  \\\n" % sam_fn
    else:
        cmd += "    -S  %s  \\\n" % sam_list_fn
    if mode == "10x":
        cmd += "    -b  %s  \\\n" % barcode_fn
    elif mode == "smartseq":
        cmd += "    -i  %s  \\\n" % barcode_fn
    else:
        cmd += "    -I  %s  \\\n" % sample
    cmd += "    -O  %s  \\\n" % out_dir
    cmd += "    -R  %s  \\\n" % snp_vcf_fn
    cmd += "    -p  %s  \\\n" % str(ncores)
    cmd += "    --minMAF  %s  \\\n" % str(min_maf)
    cmd += "    --minCOUNT  %s  \\\n" % str(min_count)
    cmd += "    --cellTAG  %s  \\\n" % str(cell_tag)
    cmd += "    --UMItag  %s  \n" % str(umi_tag)

    with open(script_fn, "w") as fp:
        fp.write(cmd)
    st = os.stat(script_fn)
    os.chmod(script_fn, st.st_mode | stat.S_IXUSR)

    # run pileup
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
    except Exception as e:
        error(str(e))
        error("Error: pileup failed (retcode '%s')." % str(ret))
        sys.exit(1)
