# baf.py - preprocess the input BAM file to generate reference-phased cell x gene AD & DP matrices.

import getopt
import gzip
import logging
import os
import pandas as pd
import random
import stat
import subprocess
import sys

from logging import error, info
from xcltk.baf.genotype import pileup
from xcltk.utils.base import assert_e, assert_n
from xcltk.utils.xlog import init_logging


def usage(fp = sys.stdout):
    s =  "\n" 
    s += "Version: %s\n" % VERSION
    s += "Usage:   %s [options]\n" % APP
    s += "\n" 
    s += "Options:\n"
    s += "  --label STR        Task label.\n"
    s += "  --sam FILE         Comma separated indexed BAM/CRAM file(s).\n"
    s += "  --samlist FILE     A file listing BAM/CRAM files, each per line.\n"
    s += "  --barcode FILE     A plain file listing all effective cell barcodes (for 10x)\n"
    s += "                     or sample IDs (for smartseq).\n"
    s += "  --snpvcf FILE      A vcf file listing all candidate SNPs.\n"
    s += "  --outdir DIR       Output dir.\n"
    s += "  --gmap FILE        Path to genetic map provided by Eagle2\n"
    s += "                     (e.g. Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz).\n"
    s += "  --eagle FILE       Path to Eagle2 binary file.\n"
    s += "  --paneldir DIR     Directory to phasing reference panel (BCF files).\n"
    s += "  --version          Print version and exit.\n"
    s += "  --help             Print this message and exit.\n"
    s += "\n"
    s += "Optional arguments:\n"
    s += "  --cellTAG STR      Cell barcode tag [%s]\n" % CELL_TAG
    s += "  --UMItag STR       UMI tag [%s]\n" % UMI_TAG
    s += "  --ncores INT       Number of threads [%d]\n" % N_CORES
    s += "  --smartseq         Run in smartseq mode.\n"
    s += "  --bulk             Run in bulk mode.\n"
    s += "\n"
    s += "Notes:\n"
    s += "  1. One and only one of `--sam` and `--samlist` should be specified.\n"
    s += "  2. For smartseq data, the order of the BAM files (in `--sam` or `--samlist`)\n"
    s += "     and the sample IDs (in `--barcode`) should match each other.\n"
    s += "\n"

    fp.write(s)


def main(argv):
    if len(argv) < 2:
        usage()
        sys.exit(1)

    init_logging(stream = sys.stdout)

    label = None
    sam_fn = sam_list_fn = barcode_fn = snp_vcf_fn = None
    out_dir = None
    gmap_fn = eagle_fn = panel_dir = None
    cell_tag, umi_tag = CELL_TAG, UMI_TAG
    ncores = N_CORES
    mode = "10x"

    opts, args = getopt.getopt(
        args = argv[1:],
        shortopts = "", 
        longopts = [
            "label=",
            "sam=", "samlist=", "barcode=", "snpvcf=",
            "outdir=",
            "gmap=", "eagle=", "paneldir=",
            "version", "help",
            
            "cellTAG=", "UMItag=",
            "ncores=",
            "smartseq", "bulk"
        ])

    for op, val in opts:
        if len(op) > 2:
            op = op.lower()
        if op in ("--label"): label = val
        elif op in ("--sam"): sam_fn = val
        elif op in ("--samlist"): sam_list_fn = val
        elif op in ("--barcode"): barcode_fn = val
        elif op in ("--snpvcf"): snp_vcf_fn = val
        elif op in ("--outdir"): out_dir = val
        elif op in ("--gmap"): gmap_fn = val
        elif op in ("--eagle"): eagle_fn = val
        elif op in ("--paneldir"): panel_dir = val
        elif op in ("--version"): error(VERSION); sys.exit(1)
        elif op in ("--help"): usage(); sys.exit(1)

        elif op in ("--celltag"): cell_tag = val
        elif op in ("--umitag"): umi_tag = val
        elif op in ("--ncores"): ncores = int(val)     # keep it in `str` format.
        elif op in ("--smartseq"): mode = "smartseq"
        elif op in ("--bulk"): mode = "bulk"
        else:
            error("invalid option: '%s'." % op)
            return(-1)

    info("xcltk BAF preprocessing starts ...")

    # check args
    info("check args ...")

    assert_n(label)
    sample = label if mode == "bulk" else None
    assert_e(gmap_fn)
    genome = "hg19" if "hg19" in gmap_fn else "hg38"
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    # other args will be checked in pileup().
    info("run in '%s' mode (genome version '%s')" % (mode, genome))

    # pileup
    info("start pileup ...")

    pileup_dir = os.path.join(out_dir, "pileup")
    if not os.path.exists(pileup_dir):
        os.mkdir(pileup_dir)
    pileup_script = os.path.join(out_dir, "run_pileup.sh")
    pileup_log_fn = os.path.join(out_dir, "pileup.log")

    pileup(
        sam_fn = sam_fn, sam_list_fn = sam_list_fn,
        barcode_fn = barcode_fn, sample = sample,
        snp_vcf_fn = snp_vcf_fn,
        out_dir = pileup_dir,
        mode = mode,
        cell_tag = cell_tag, umi_tag = umi_tag,
        ncores = ncores,
        min_maf = 0.1, min_count = 20,
        script_fn = pileup_script,
        log_fn = pileup_log_fn
    )


APP = "baf.py"
VERSION = "0.0.1"

CELL_TAG = "CB"
N_CORES = 1
UMI_TAG = "UB"


if __name__ == "__main__":
    main(sys.argv)