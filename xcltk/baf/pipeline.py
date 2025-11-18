# pipeline.py - preprocess the input BAM file to generate reference-phased cell x gene AD & DP matrices.



import getopt
import os
import sys

from logging import error, info
from .fc.main import afc_wrapper as baf_fc      # BAF feature counting
from .genotype import pileup, vcf_add_genotype
from .refphase import ref_phasing
from ..config import APP, VERSION
from ..utils.base import assert_e, assert_n
from ..utils.vcf import vcf_index, vcf_merge, vcf_split_chrom
from ..utils.xlog import init_logging



COMMAND = "baf"

CELL_TAG = "CB"
UMI_TAG = "UB"
MIN_COUNT = 11
MIN_MAF = 0.1
N_CORES = 1



def usage(fp = sys.stdout):
    s =  "\n" 
    s += "Version: %s\n" % VERSION
    s += "Usage:   %s %s [options]\n" % (APP, COMMAND)
    s += "\n" 
    s += "Options:\n"
    s += "  --label STR        Task label.\n"
    s += "  --sam FILE         Comma separated indexed BAM/CRAM file(s).\n"
    s += "  --samList FILE     A file listing BAM/CRAM files, each per line.\n"
    s += "  --barcode FILE     A plain file listing all effective cell barcodes, for\n"
    s += "                     droplet-based data, e.g., 10x Genomics.\n"
    s += "  --sampleList FILE  A plain file listing sample IDs, one ID per BAM, for\n"
    s += "                     well-based data, e.g., SMART-seq.\n"
    s += "  --snpvcf FILE      A VCF file listing all candidate SNPs.\n"
    s += "  --region FILE      A TSV file listing target features. The first 4 columns are:\n"
    s += "                     chrom, start, end (both 1-based and inclusive), name.\n"
    s += "  --outdir DIR       Output dir.\n"
    s += "  --gmap FILE        Path to genetic map provided by Eagle2\n"
    s += "                     (e.g. Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz).\n"
    s += "  --eagle FILE       Path to Eagle2 binary file.\n"
    s += "  --paneldir DIR     Directory to phasing reference panel (BCF files).\n"
    s += "  --version          Print version and exit.\n"
    s += "  --help             Print this message and exit.\n"
    s += "\n"
    s += "Optional arguments:\n"
    s += "  --refCell FILE     A plain file listing reference cells, one per line.\n"
    s += "  --cellTAG STR      Cell barcode tag; Set to None if not available [%s]\n" % CELL_TAG
    s += "  --UMItag STR       UMI tag; Set to None if not available [%s]\n" % UMI_TAG
    s += "  --minCOUNT INT     Mininum aggragated count for SNP [%d]\n" % MIN_COUNT
    s += "  --minMAF FLOAT     Mininum minor allele fraction for SNP [%f]\n" % MIN_MAF
    s += "  --ncores INT       Number of threads [%d]\n" % N_CORES
    s += "\n"
    s += "Notes:\n"
    s += "1. One and only one of `--sam` and `--samlist` should be specified.\n"
    s += "2. For well-based data, the order of the BAM files (in `--sam` or `--samlist`)\n"
    s += "   and the sample IDs (in `--sampleList`) should match each other.\n"
    s += "3. For bulk data, the label (`--label`) will be used as the sample ID.\n"
    s += "\n"

    fp.write(s)


    
def pipeline_main(argv):
    if len(argv) <= 2:
        usage()
        sys.exit(0)

    init_logging(stream = sys.stdout)

    label = None
    sam_fn = sam_list_fn = None
    barcode_fn = sample_id_fn = None
    snp_vcf_fn = region_fn = None
    out_dir = None
    gmap_fn = eagle_fn = panel_dir = None
    ref_cell_fn = None
    cell_tag, umi_tag = CELL_TAG, UMI_TAG
    min_count, min_maf = MIN_COUNT, MIN_MAF
    ncores = N_CORES

    opts, args = getopt.getopt(
        args = argv[2:],
        shortopts = "", 
        longopts = [
            "label=",
            "sam=", "samList=", 
            "barcode=", "sampleList=",
            "snpvcf=", "region=",
            "outdir=",
            "gmap=", "eagle=", "paneldir=",
            "version", "help",

            "refCell=",
            "cellTAG=", "UMItag=",
            "minCOUNT=", "minMAF=",
            "ncores="
        ])

    for op, val in opts:
        if len(op) > 2:
            op = op.lower()
        if op in ("--label"): label = val
        elif op in ("--sam"): sam_fn = val
        elif op in ("--samlist"): sam_list_fn = val
        elif op in ("--barcode"): barcode_fn = val
        elif op in ("--samplelist"): sample_id_fn = val
        elif op in ("--snpvcf"): snp_vcf_fn = val
        elif op in ("--region"): region_fn = val
        elif op in ("--outdir"): out_dir = val
        elif op in ("--gmap"): gmap_fn = val
        elif op in ("--eagle"): eagle_fn = val
        elif op in ("--paneldir"): panel_dir = val
        elif op in ("--version"): sys.stdout.write(VERSION + "\n"); sys.exit(0)
        elif op in ("--help"): usage(); sys.exit(0)

        elif op in ("--refcell"): ref_cell_fn = val
        elif op in ("--celltag"): cell_tag = val
        elif op in ("--umitag"): umi_tag = val
        elif op in ("--mincount"): min_count = int(val)
        elif op in ("--minmaf"): min_maf = float(val)
        elif op in ("--ncores"): ncores = int(val)     # keep it in `str` format.
        else:
            error("invalid option: '%s'." % op)
            return(-1)
        
    ret = pipeline_wrapper(
        label = label,
        sam_fn = sam_fn, sam_list_fn = sam_list_fn, 
        barcode_fn = barcode_fn, sample_id_fn = sample_id_fn,
        snp_vcf_fn = snp_vcf_fn, region_fn = region_fn,
        out_dir = out_dir,
        gmap_fn = gmap_fn, eagle_fn = eagle_fn, panel_dir = panel_dir,
        ref_cell_fn = ref_cell_fn,
        cell_tag = cell_tag, umi_tag = umi_tag,
        min_count = min_count, min_maf = min_maf,
        ncores = ncores
    )
    
    info("All Done!")

    return(ret)
        

    
def pipeline_wrapper(
    label,
    sam_fn = None, sam_list_fn = None, 
    barcode_fn = None, sample_id_fn = None,
    snp_vcf_fn = None, region_fn = None,
    out_dir = None,
    gmap_fn = None, eagle_fn = None, panel_dir = None,
    ref_cell_fn = None,
    cell_tag = "CB", umi_tag = "UB",
    min_count = 11, min_maf = 0.1,
    ncores = 1
):
    info("xcltk BAF preprocessing starts ...")

    # check args
    info("check args ...")

    mode = None
    sample_id = None

    assert_n(label)
    assert (not sam_fn) ^ (not sam_list_fn)
    assert not (barcode_fn and sample_id_fn)
    if barcode_fn:
        assert_e(barcode_fn)
        mode = "droplet"
    elif sample_id_fn:
        assert_e(sample_id_fn)
        mode = "well"
    else:
        mode = "bulk"
        sample_id = label

    assert_e(snp_vcf_fn)
    assert_e(region_fn)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok = True)
    script_dir = os.path.join(out_dir, "scripts")
    os.makedirs(script_dir, exist_ok = True)

    assert_e(gmap_fn)
    assert_e(eagle_fn)
    assert_e(panel_dir)
    for chrom in range(1, 23):
        assert_e(os.path.join(panel_dir, "chr%d.genotypes.bcf" % chrom))
        assert_e(os.path.join(panel_dir, "chr%d.genotypes.bcf.csi" % chrom))
        
    if ref_cell_fn is not None:
        assert_e(ref_cell_fn)

    genome = "hg19" if "hg19" in gmap_fn else "hg38"

    # other args will be checked in pileup() and ref_phasing().
    info("run in '%s' mode (genome version '%s')" % (mode, genome))

    step = 1


    # pileup
    info("start pileup ...")

    pileup_dir = os.path.join(out_dir, "%d_pileup" % step)
    if not os.path.exists(pileup_dir):
        os.makedirs(pileup_dir, exist_ok = True)
    pileup_script_dir = os.path.join(script_dir, "pileup")
    os.makedirs(pileup_script_dir, exist_ok = True)
    pileup_script = os.path.join(pileup_script_dir, "run_pileup.sh")
    pileup_log_fn = os.path.join(pileup_script_dir, "pileup.log")

    pileup_vcf_fn, p_raw, p_new = pileup(
        sam_fn = sam_fn, sam_list_fn = sam_list_fn,
        barcode_fn = barcode_fn, sample_id_fn = sample_id_fn, 
        sample_id = sample_id,
        snp_vcf_fn = snp_vcf_fn,
        out_dir = pileup_dir,
        mode = mode,
        cell_tag = cell_tag, umi_tag = umi_tag,
        ncores = ncores,
        min_count = min_count, min_maf = min_maf,
        script_fn = pileup_script,
        log_fn = pileup_log_fn
    )
    
    info("pileup #SNP raw=%d; post-filtering=%d." % (p_raw, p_new))

    assert_e(pileup_vcf_fn)

    info("pileup VCF is '%s'." % pileup_vcf_fn)
    step += 1


    # prepare VCF files for phasing
    info("prepare VCF files for phasing ...")

    phasing_dir = os.path.join(out_dir, "%d_phasing" % step)
    if not os.path.exists(phasing_dir):
        os.makedirs(phasing_dir, exist_ok = True)
    phasing_chrom_dir = os.path.join(phasing_dir, "chrom_phasing")
    os.makedirs(phasing_chrom_dir, exist_ok = True)
    phasing_script_dir = os.path.join(script_dir, "phasing")
    os.makedirs(phasing_script_dir, exist_ok = True)
    phasing_script_prefix = os.path.join(phasing_script_dir, "run_phasing")
    phasing_log_prefix = os.path.join(phasing_script_dir, "phasing")

    
    # add genotypes
    info("add genotypes ...")

    genotype_vcf_fn = os.path.join(phasing_dir, "%s.genotype.vcf.gz" % label)
    vcf_add_genotype(
        in_fn = pileup_vcf_fn,
        out_fn = genotype_vcf_fn,
        sample = label,
        chr_prefix = True,     # add "chr" prefix
        sort = True,
        unique = True
    )

    
    # split VCF by chromosomes.
    info("split VCF by chromosomes ...")

    valid_chroms = []       # has "chr" prefix
    target_vcf_list = []
    res = vcf_split_chrom(
        fn = genotype_vcf_fn,
        out_dir = phasing_chrom_dir,
        label = label,
        chrom_list = ["chr" + str(i) for i in range(1, 23)],
        out_prefix_list = None,
        verbose = True
    )
    for chrom, n_variants, vcf_fn in res:
        if n_variants > 0:
            valid_chroms.append(chrom)
            target_vcf_list.append(vcf_fn)
            vcf_index(vcf_fn)

    info("%d chromosome VCFs are outputted with variants." % len(valid_chroms))
    #os.remove(genotype_vcf_fn)

    
    # reference phasing
    info("reference phasing ...")

    ref_vcf_list = [os.path.join(panel_dir, "%s.genotypes.bcf" % chrom) \
                    for chrom in valid_chroms]
    out_prefix_list = [os.path.join(phasing_chrom_dir, "%s_%s.phased" % \
                    (label, chrom)) for chrom in valid_chroms]
    ref_phasing(
        target_vcf_list = target_vcf_list,
        ref_vcf_list = ref_vcf_list,
        out_prefix_list = out_prefix_list,
        gmap_fn = gmap_fn,
        eagle_fn = eagle_fn,
        out_dir = phasing_chrom_dir,
        ncores = ncores,
        script_fn_prefix = phasing_script_prefix,
        log_fn_prefix = phasing_log_prefix,
        verbose = True
    )

    info("phased VCFs are in dir '%s'." % phasing_chrom_dir)

    
    # merge phased VCFs
    info("merge phased VCFs ...")

    phased_vcf_list = ["%s.vcf.gz" % prefix for prefix in out_prefix_list]
    for fn in phased_vcf_list:
        assert_e(fn)

    phased_vcf_fn = os.path.join(phasing_dir, "%s.phased.vcf.gz" % label)
    vcf_merge(phased_vcf_list, phased_vcf_fn, sort = True)

    info("merged VCF is '%s'." % phased_vcf_fn)
    step += 1


    # allele-specific feature counting.
    info("BAF feature counting ...")

    fc_dir = os.path.join(out_dir, "%d_baf_fc" % step)
    if not os.path.exists(fc_dir):
        os.makedirs(fc_dir, exist_ok = True)

    baf_fc(
        sam_fn = sam_fn, 
        barcode_fn = barcode_fn,
        region_fn = region_fn, 
        phased_snp_fn = phased_vcf_fn,
        out_dir = fc_dir,
        sam_list_fn = sam_list_fn,
        sample_ids = sample_id,
        sample_id_fn = sample_id_fn,
        debug_level = 0,
        ncores = ncores,
        cellsnp_dir = pileup_dir, 
        ref_cell_fn = ref_cell_fn,
        cell_tag = cell_tag, umi_tag = umi_tag,
        min_count = 1, min_maf = 0,
        output_all_reg = True, no_dup_hap = True,
        min_mapq = 20, min_len = 30,
        incl_flag = 0, excl_flag = None,
        no_orphan = True
    )

    info("feature BAFs are at '%s'." % fc_dir)
    step += 1
