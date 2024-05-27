# vcf.py - VCF file routine

import gzip
import numpy as np
import os
import pandas as pd
import random
import subprocess
import sys
import time

from logging import error, info
from logging import warning as warn
from .base import assert_e, assert_n


def vcf_load(fn):
    """Load VCF file
    
    Load the header and variants part of the VCF file.

    Parameters
    ----------
    fn : str
        The VCF file.

    Returns
    -------
    Tuple
        A tuple of two elements:
        (1) variants : pandas DataFrame object;
            The variants part of the VCF. 
            Its column names are determined by the last line of the 
            `header` part.
        (2) header : list;
            The header part of the VCF.
            One line (without tailing line separator) per element in the list.
    """
    # load header
    fp = None
    if fn.lower().endswith(".gz"):
        fp = gzip.open(fn, "rt")
    else:
        fp = open(fn, "r")
        
    header = []
    pre_line = None
    for line in fp:
        if not line or line[0] != "#":
            break
        pre_line = line
        header.append(line.rstrip())
    
    fp.close()
    
    if not pre_line:
        raise IOError()
    assert len(pre_line) > 6
    assert pre_line[:6] == "#CHROM"
    
    # load variants
    variants = pd.read_csv(fn, sep = "\t", header = None, comment = "#",
                            dtype = {0: str})
    variants.columns = pre_line.strip()[1:].split("\t")
    
    return((variants, header))


def vcf_save(variants, header, fn, is_gzip = None):
    """Save VCF file

    Save the header and variants into VCF file, in either plain VCF or 
    BGZF compressed format.
    When saving in BGZF format, it will use the `pysam` python package (if 
    available) or call the `bgzip` command-line tool (otherwise).

    Parameters
    ----------
    variants : pandas DataFrame
        The variants part of the VCF file. A pandas DataFrame object inherited
        from the one returned by `vcf_load()`.
    header : list
        The header part of the VCF file.
        One VCF header line (without tailing line separator) per element in
        the list.
    fn : str
        The output VCF file.
    is_gzip : bool
        Whether to output in BGZF format. If `None`, the value will be 
        determined by checking the file name suffix of `fn`. 
        Specifically, `True` if `fn` ends with ".gz" or ".GZ", 
        otherwise `False`.

    Returns
    -------
    None
    """
    if is_gzip is None:
        is_gzip = fn.lower().endswith(".gz")

    assert len(header) > 0
    assert header[-1].startswith("#CHROM")
    header_str = "\n".join(header) + "\n"

    # note pandas `to_csv` will always output in bytes/binary format
    # if the file has suffix ".gz".
    tmp_fn = __gen_tmp_filename(fn, "vcf_save")
    with open(tmp_fn, "w") as fp:
        fp.write(header_str)
    variants.to_csv(
        tmp_fn, sep = "\t", header = False, index = False, mode = "a")
    
    if is_gzip:
        fn_zip = vcf_bgzip(tmp_fn, fn, is_in_gzip = False)
        if fn_zip != fn:
            os.rename(fn_zip, fn)
        os.remove(tmp_fn)
    else:
        os.rename(tmp_fn, fn)


def vcf_bgzip(in_fn, out_fn = None, is_in_gzip = False):
    assert_e(in_fn)
    if not out_fn:
        if in_fn.lower().endswith(".gz"):
            out_fn = in_fn[:-3] + ".vcf_bgzip.gz"
            warn("input file has suffix '.gz'; write to '%s'" % out_fn)
        else:
            out_fn = in_fn + ".gz"
    try:
        import pysam
    except:
        try:
            proc = subprocess.Popen(
                args = "bgzip -c %s > %s" % (in_fn, out_fn),
                shell = True,
                executable = "/bin/bash",
                stdout = subprocess.PIPE,
                stderr = subprocess.PIPE
            )
            outs, errs = proc.communicate()
            ret = proc.returncode
            if ret != 0:
                raise RuntimeError(str(errs.decode()))
        except:
            error("Error: bgzip failed (retcode '%s')." % str(ret))
            sys.exit(1)        
    else:
        in_fp = None
        if is_in_gzip:
            in_fp = gzip.open(in_fn, "rt")
        else:
            in_fp = open(in_fn, "r")
        out_fp = pysam.BGZFile(out_fn, "w")
        for line in in_fp:
            out_fp.write(line.encode())
        in_fp.close()
        out_fp.close()
    return(out_fn)


def vcf_index(fn, idx_fn = None, idx_type = None, ncores = 1):
    assert_e(fn)
    if idx_type is None:
        if idx_fn and idx_fn.endswith(".tbi"):
            idx_type = "tbi"
        else:
            idx_type = "csi"
    assert idx_type in ("csi", "tbi")

    if idx_fn is None:
        idx_fn = "%s.%s" % (fn, idx_type)

    try:
        import pysam
    except:
        cmd = None
        if idx_type == "csi":
            cmd = "bcftools index -c -o %s --threads %d %s" % \
                (idx_fn, ncores, fn)
        else:
            cmd = "bcftools index -t -o %s --threads %d %s" % \
                (idx_fn, ncores, fn)
        try:
            proc = subprocess.Popen(
                args = cmd,
                shell = True, 
                executable = "/bin/bash", 
                stdout = subprocess.PIPE, 
                stderr = subprocess.PIPE
            )
            outs, errs = proc.communicate()
            ret = proc.returncode
            if ret != 0:
                raise RuntimeError(str(errs.decode()))
        except:
            error("Error: bcftools index failed (retcode '%s')." % str(ret))
            sys.exit(1)
    else:
        csi = True if idx_type == "csi" else False
        pysam.tabix_index(
            fn, preset = "vcf", index = idx_fn, 
            force = True, csi = csi)
    return(idx_fn)


def vcf_merge(in_fn_list, out_fn, sort = False):
    assert len(in_fn_list) > 0
    header = None
    all_vars = None
    for i, fn in enumerate(in_fn_list):
        if i == 0:
            variants, header = vcf_load(fn)
            all_vars = variants
        else:
            variants, _ = vcf_load(fn)
            all_vars = pd.concat([all_vars, variants], ignore_index = True)
    if sort:
        all_vars = all_vars.sort_values(by = ["CHROM", "POS", "REF", "ALT"])
    vcf_save(all_vars, header, out_fn)
    

# for processing the "chr" prefix with command-line tools, please
# refer to:
# 1. bcftools reheader; https://www.biostars.org/p/198660/
# 2. bcftools index -s && bcftools annotate --rename-chrs;
#    https://github.com/samtools/bcftools/issues/747
# 3. bcftools; https://alkesgroup.broadinstitute.org/Eagle/

def vcf_add_chr_prefix(in_fn, out_fn):
    variants, header = vcf_load(in_fn)
    variants, header = vcf_add_chr_prefix_core(variants, header)
    vcf_save(variants, header, out_fn)


def vcf_add_chr_prefix_core(variants, header, inplace = False):
    if not inplace:
        variants = variants.copy()
        header = header.copy()
    assert "CHROM" in variants.columns
    if np.any(~(variants["CHROM"].astype(str).str.startswith("chr"))):
        format_chrom = lambda x: x if x.startswith("chr") else "chr" + x
        variants["CHROM"] = variants["CHROM"].astype(str).map(format_chrom)
    for i in range(len(header)):
        h = header[i]
        if h.startswith("##contig=") and not h.startswith("##contig=<ID=chr"):
            h = h.replace("##contig=<ID=", "##contig=<ID=chr")
        header[i] = h
    return(variants, header)


def vcf_remove_chr_prefix(in_fn, out_fn):
    variants, header = vcf_load(in_fn)
    variants, header = vcf_remove_chr_prefix_core(variants, header)
    vcf_save(variants, header, out_fn)


def vcf_remove_chr_prefix_core(variants, header, inplace = False):
    if not inplace:
        variants = variants.copy()
        header = header.copy()
    assert "CHROM" in variants.columns
    if np.any(variants["CHROM"].astype(str).str.startswith("chr")):
        format_chrom = lambda x: x[3:] if x.startswith("chr") else x
        variants["CHROM"] = variants["CHROM"].astype(str).map(format_chrom)
    header = [s.replace("##contig=<ID=chr", "##contig=<ID=") for s in header]
    return(variants, header)


def vcf_hdr_check_contig(in_fn, out_fn = None):
    assert_e(in_fn)
    if out_fn is None:
        out_fn = in_fn
    variants, header = vcf_load(in_fn)
    variants, header = vcf_hdr_check_contig_core(variants, header)
    vcf_save(variants, header, out_fn)


def vcf_hdr_check_contig_core(variants, header, inplace = False):
    if not inplace:
        variants = variants.copy()
        header = header.copy()
    assert "CHROM" in variants.columns
    chrom_list = np.sort(variants["CHROM"].unique())
    assert len(chrom_list) > 0

    hdr_contig = []
    hdr_oth = []
    idx_first_contig = None
    chrom_exists = [False] * len(chrom_list)
    for idx_line, line in enumerate(header):
        if line.startswith("##contig=<ID="):
            hdr_contig.append(line)
            if idx_first_contig is None:
                idx_first_contig = idx_line
        else:
            hdr_oth.append(line)
        for idx_chrom, chrom in enumerate(chrom_list):
            if line.startswith("##contig=<ID=%s" % chrom):
                chrom_exists[idx_chrom] = True

    for exist, chrom in zip(chrom_exists, chrom_list):
        if not exist:
            hdr_contig.append("##contig=<ID=%s>" % chrom)

    if idx_first_contig is None:
        header = hdr_oth[:-1] + sorted(hdr_contig) + [hdr_oth[-1]]
    else:
        header = hdr_oth[:idx_first_contig] + sorted(hdr_contig) + \
            hdr_oth[idx_first_contig:]

    return(variants, header)
    

def vcf_split_chrom(
    fn, out_dir, 
    label = None, 
    chrom_list = None,
    out_prefix_list = None,
    verbose = False
):
    """Split VCF file by chromosomes

    Split VCF file by chromosomes, one output VCF per chromosome.

    Parameters
    ----------
    fn : str
        The input VCF file.
    out_dir : str
        The output dir to store the output chromosome-specific VCF files.
    label : str
        The task label. It not `None`, it will be part of the prefix of output
        files (e.g., <out_dir>/<label>_chr<chrom>.vcf.gz).
        Otherwise, the `output_prefix_list` will be used (e.g., 
        <out_dir/<out_prefix_list[0]>.vcf.gz).
    chrom_list : list
        The list of chromosomes in each output VCF file.
        If `None`, it will use all unique chromosomes in the input VCF file.
    out_prefix_list : list
        The list of prefix of the output VCF files. If not `None`, its length
        should be the same with `chrom_list`. See `label` for details.
    verbose: bool
        Whether to show detailed logging information.
    
    Returns
    -------
    list
        A list of tuples. Each tuple has three elements: 
        (1) chrom : str;
            Chromosome name.
        (2) n_variants : int;
            Number of variants on the chromosome.
        (3) vcf_fn : str;
            The output VCF file.

    Raises
    ------
    AssertionError
        When some arguments are invalid.
    """
    assert_e(fn)
    assert_n(out_dir)
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    variants, header = vcf_load(fn)
    assert variants is not None and variants.shape[0] > 0
    assert len(header) > 0
    if chrom_list is None:
        chrom_list = variants["CHROM"].unique()

    if label is None:
        assert out_prefix_list is not None and \
            len(out_prefix_list) == len(chrom_list)
    else:
        assert out_prefix_list is None
        format_chrom = lambda x: x if x.startswith("chr") else "chr" + x
        out_prefix_list = ["%s_%s" % (label, format_chrom(chrom)) \
            for chrom in chrom_list]
    
    if verbose:
        info("split %d chromosomes ..." % len(chrom_list))

    out_fn_list = []
    for chrom, out_prefix in zip(chrom_list, out_prefix_list):
        if verbose:
            info("process chrom '%s' ..." % chrom)
        variants_chrom = variants.loc[variants["CHROM"] == chrom, :]
        n_variants = variants_chrom.shape[0]
        if n_variants <= 0:
            warn("no any variants for chrom '%s'." % chrom)
        out_fn = os.path.join(out_dir, "%s.vcf.gz" % out_prefix)
        vcf_save(variants_chrom, header, out_fn)
        out_fn_list.append((chrom, n_variants, out_fn))

    return(out_fn_list)


def __gen_tmp_filename(fn, func = None, suffix = None):
    items = [fn]
    if func:
        items.append(func)
    items.append(time.strftime("%Y%m%d", time.localtime()))
    items.append(str(random.randint(1000, 10000)))
    if suffix:
        items.append(suffix)
    tmp_fn = "_".join(items) + ".tmp"
    return(tmp_fn)