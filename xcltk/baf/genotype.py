# genotype.py - genotyping SNPs.


import os
import stat
import subprocess
import sys

from logging import error, info
from logging import warning as warn
from ..utils.base import assert_e, assert_n
from ..utils.csp_io import load_data as csp_load_data
from ..utils.csp_io import save_data as csp_save_data
from ..utils.vcf import vcf_load, vcf_save, \
    vcf_add_chr_prefix_core, vcf_remove_chr_prefix_core, \
    vcf_hdr_check_contig_core



def pileup(
    sam_fn = None, sam_list_fn = None, 
    barcode_fn = None, sample_id_fn = None,
    sample_id = None,
    snp_vcf_fn = None,
    out_dir = None,
    mode = "droplet",
    cell_tag = "CB", umi_tag = "UB",
    ncores = 1,
    min_count = 11, min_maf = 0.1,
    script_fn = None, log_fn = None
):
    """Pileup indexed BAM file, supporting both single-cell and bulk data.

    The function takes as input the BAM file and a list of SNPs for pileup,
    and outputs three sparse matrices, "AD" (alternative allele counts), 
    "DP" (reference and alternative allele counts), and "OTH" 
    (counts of other alleles).
    For pileup, it internally calls the `cellsnp-lite` command-line tool.

    Parameters
    ----------
    sam_fn : str
        Comma separated indexed BAM/CRAM file(s).
    sam_list_fn : str
        A file listing BAM/CRAM files, each per line.
        Note that The input BAM file(s) should be specified by one and only
        one of `sam_fn` and `sam_list_fn`.
    barcode_fn : str
        A plain file listing all effective cell barcodes (for droplet-based
        data, e.g., 10x Genomics).
    sample_id_fn : str
        A plain file listing all sample IDs (for well-based data, e.g.,
        SMART-seq).
        For well-based data, the order of the BAM files (in `sam_fn` or 
        `sam_list_fn`) and the sample IDs should match each other.
    sample_id : str
        Sample ID (for bulk data).
    snp_vcf_fn : str
        A vcf file listing all candidate SNPs.
    out_dir : str
        Output dir.
    mode : str
        One of "droplet", "well", and `bulk`.
    cell_tag : str
        Cell barcode tag, set to `None` to turn it off.
    umi_tag : str
        UMI tag, set to `None` to turn it off.
    ncores : int
        Number of threads.
    min_count : int
        Minimum aggregated count.
    min_maf : float
        Minimum minor allele frequency.
    script_fn : str
        Path to the script file that runs cellsnp-lite. If `None`, use default
        path `<out_dir>/run_pileup.sh`.
    log_fn : str
        Path to the logging file that records the output of cellsnp-lite. If
        `None`, use default path `<out_dir>/pileup.log`.

    Returns
    -------
    str
        The VCF file containing final output SNPs.

    Raises
    ------
    AssertionError
        when some input file or dirs are invalid.
    """
    # check args
    assert mode in ("droplet", "well", "bulk")

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


    if mode == "droplet":
        assert_e(barcode_fn)
    elif mode == "well":
        assert_e(sample_id_fn)
        sid_list = [line.strip() for line in open(sample_id_fn, "r")]
        assert len(sam_list) == len(sid_list)
    else:
        assert_n(sample_id)
        assert len(sam_list) == 1

    assert_e(snp_vcf_fn)

    if mode in ("well", "bulk"):
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
        os.makedirs(out_dir, exist_ok = True)
    raw_dir = os.path.join(out_dir, 'raw')
    if not os.path.exists(raw_dir):
        os.makedirs(raw_dir, exist_ok = True)
    
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
    if mode == "droplet":
        cmd += "    -b  %s  \\\n" % barcode_fn
    elif mode == "well":
        cmd += "    -i  %s  \\\n" % sample_id_fn
    else:
        cmd += "    -I  %s  \\\n" % sample_id
    cmd += "    -O  %s  \\\n" % raw_dir
    cmd += "    -R  %s  \\\n" % snp_vcf_fn
    cmd += "    -p  %s  \\\n" % str(ncores)
    cmd += "    --minCOUNT  %s  \\\n" % str(min_count)
    cmd += "    --minMAF  %s  \\\n" % str(min_maf)
    cmd += "    --cellTAG  %s  \\\n" % str(cell_tag)
    cmd += "    --UMItag  %s  \\\n" % str(umi_tag)
    cmd += "    --gzip    \n"

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
        if ret != 0:
            raise RuntimeError(str(errs.decode()))
    except Exception as e:
        error(str(e))
        error("Error: pileup failed (retcode '%s')." % str(ret))
        sys.exit(1)


    out_vcf_fn, p_raw, p_new = filter_snps(
        in_dir = raw_dir,
        out_dir = out_dir,
        min_count = min_count,
        min_maf = min_maf
    )
    return(out_vcf_fn, p_raw, p_new)

        

def filter_snps(in_dir, out_dir, min_count, min_maf):
    """Further filter SNPs.
    
    Filter SNPs given `minMAF` and `minCOUNT`, considering only REF and ALT
    (AD & DP) alleles, but not OTH alleles.
    Cellsnp-lite (at least v1.2.3 and before) may keep some SNPs unexpectly,
    e.g., SNP with `AD=9;DP=9;OTH=1` when `minMAF=0.1; minCOUNT=10`.
    """
    adata = csp_load_data(in_dir)
    n, p = adata.shape
    adata.var.index = [str(i) for i in range(p)]
    
    df = adata.var.copy()
    df["AD"] = df["INFO"].map(
        lambda x: int(x.split(";")[0].split("=")[1])
    )
    df["DP"] = df["INFO"].map(
        lambda x: int(x.split(";")[1].split("=")[1])
    )
    df = df[df["DP"] >= min_count].copy()
    df["BAF"] = df["AD"] / df["DP"]
    df = df[(df["BAF"] >= min_maf) & (df["BAF"] <= 1-min_maf)].copy()
    
    adata = adata[:, df.index].copy()
    n_new, p_new = adata.shape
    
    # save SNPs.
    csp_save_data(adata, out_dir)
    out_vcf_fn = os.path.join(out_dir, "cellSNP.base.vcf.gz")
    return(out_vcf_fn, p, p_new)



def vcf_add_genotype(
    in_fn, out_fn, 
    sample, 
    chr_prefix = None, 
    sort = True,
    unique = True
):
    """Add genotypes (`GT` field) in VCF file

    The function adds variant genotpyes (i.e., the `GT` field) into the 
    VCF file.
    Currently, all genotypes will be simply set as "0/1" assuming all input
    variants are heterzygous.

    Parameters
    ----------
    in_fn : str
        The input VCF file.
    out_fn : str
        The oupput VCF file.
    sample : str
        The sample name, used in the sample field of the output VCF.
    chr_prefix : bool
        Should the chromosome names have "chr" prefix. 
        If `None`, keep it unchanged.
    sort : bool
        Whether the variants should be sorted. If `True`, the SNPs will be
        sorted based on `CHROM`, `POS`, `REF`, `ALT`.
    unique: bool
        Whether the variants should be unique. 
        If `True`, the duplicate variants (based on `CHROM`, `POS`, `REF`,
        `ALT`) will be discarded, keeping the first unique records only.

    Returns
    -------
    None

    Raises
    ------
    ValueError
    AssertionError
    """
    variants, header = vcf_load(in_fn)

    if variants is None or variants.shape[0] <= 0:
        raise ValueError("no records in vcf '%s'." % in_fn)
    assert len(variants.columns) >= 8
    assert variants.columns[0] == "CHROM"

    
    # add genotype
    if "FORMAT" in variants.columns:
        warn("FORMAT in vcf '%s'." % in_fn)
    if sample in variants.columns:
        raise ValueError("sample name '%s' in vcf '%s'." % (sample, in_fn))

    header[-1] = header[-1] + "\tFORMAT\t%s" % sample
    variants["FORMAT"] = "GT"
    variants[sample] = "0/1"

    
    # make sure the vcf header "contig" lines are complete.
    variants, header = vcf_hdr_check_contig_core(variants, header)

    
    # process "chr" prefix
    if chr_prefix is not None:
        if chr_prefix:
            variants, header = vcf_add_chr_prefix_core(variants, header)
        else:
            variants, header = vcf_remove_chr_prefix_core(variants, header)

            
    # sort SNPs and drop duplicates
    if sort:
        variants = variants.sort_values(by = ["CHROM", "POS", "REF", "ALT"])
    if unique:
        variants = variants.drop_duplicates(
            subset = ["CHROM", "POS", "REF", "ALT"])

    vcf_save(variants, header, out_fn)
