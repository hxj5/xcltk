# fixref.py - fix reference allele mismatches of SNPs; mimic `bcftools +fixref`.

#this script is aimed to make sure REFs match certain genome reference build.
#it would change corresponding ALT & GT while delete other fields in FORMAT

#features:
#  support gzip/bgzip input and output vcf
#  support multiple alt alleles
#  support multiple samples
#  support only ploidy = 2 now
#  indels would be filtered

#TODO:
# - support other ploidy
# - implement fixref_wrapper().


import getopt
import pysam
import sys

from ..config import APP, VERSION

COMMAND = "fixref"


def __fix_rec(ref0, rec, out, verbose = False):
    """Check and fix reference for one VCF record.
    
    Check and fix reference for one VCF record, updating corresponding 
    ALT and GT.

    Parameters
    ----------
    ref0 : str
        Real REF from fasta.
    rec : pysam.VariantRecord
        The VCF record whose REF is to be checked.
    out : pysam.VariantFile
        The output VCF object.
    verbose : bool
        Whether to print detailed log info.

    Returns
    -------
    int
        The running state:
        * -1, if error;
        * 0, if ref == ref0, i.e., no fix is needed
        * 1, if ref != ref0 and the VCF record is successfully fixref-ed.
    pysam.VariantRecord
        The checked VCF record.
    """
    for a in rec.alleles:
        if len(a) != 1 or a not in "ACGTN":
            if verbose:
                sys.stderr.write("[W::%s] %s:%d skipped for invalid alleles\n" % (COMMAND, rec.chrom, rec.pos))
            return((-1, None))

    if "GT" not in rec.format: 
        if verbose:
            sys.stderr.write("[W::%s] %s:%d skipped for missing GT\n" % (COMMAND, rec.chrom, rec.pos))
        return((-1, None))

    new_alleles = []
    new_samples = []

    if rec.ref == ref0:
        new_alleles = rec.alleles
        for sample in rec.samples.values():
            sep = "|" if sample.phased else "/"
            ra_idx = sample.get("GT", ())      # index of ref/alt, a tuple
            gt = sep.join([str(a) for a in ra_idx if a is not None])
            if not gt: gt = "."
            new_samples.append({
                "OGT":gt,
                "GT":ra_idx
            })
    else:
        alleles = set()
        sample_alleles = []
        old_gt = []
        sep = None
        for i, sample in enumerate(rec.samples.values()):
            ra_idx = sample.get("GT", ())
            sep = "|" if sample.phased else "/"
            gt = sep.join([str(a) for a in ra_idx if a is not None])
            if not gt: gt = "."
            old_gt.append(gt)
            if len(ra_idx) == 2:
                idx1, idx2 = ra_idx[:2]
                try:
                    allele1 = rec.alleles[idx1]
                    allele2 = rec.alleles[idx2]
                except IndexError:
                    if verbose:
                        sys.stderr.write("[W::%s] %s:%d skipped for invalid alleles\n" % (COMMAND, rec.chrom, rec.pos))
                    return((-1, None))
                alleles.update([allele1, allele2])
                sample_alleles.append([allele1, allele2])
            elif len(ra_idx) == 1 and ra_idx[0] is None:      # GT = "."
                sample_alleles.append([])
            else:
                if verbose:
                    sys.stderr.write("[W::%s] %s:%d skipped for invalid GT in Sample %d\n" % (COMMAND,
                                     rec.chrom, rec.pos, i + 1))
                return((-1, None))

        alleles.discard(ref0)
        new_alts = [a for a in rec.alleles if a in alleles]   # try to keep original order of alleles
        if len(new_alts) != len(alleles):
            if verbose:
                sys.stderr.write("[W::%s] %s:%d skipped for invalid alleles\n" % (COMMAND, rec.chrom, rec.pos))
            return((-1, None))

        if sep is None:
            if verbose:
                sys.stderr.write("[W::%s] %s:%d skipped for none sample\n" % (COMMAND, rec.chrom, rec.pos))
            return((-1, None))

        new_alleles = [ref0] + new_alts
        for i, smp_ale in enumerate(sample_alleles):
            new_gt = ()
            if len(smp_ale) > 0:
                idx1, idx2 = new_alleles.index(smp_ale[0]), new_alleles.index(smp_ale[1])
                new_gt = (idx1, idx2)
            new_samples.append({
                "OGT":old_gt[i],
                "GT":new_gt
            })

    ret = 0 if rec.ref == ref0 else 1

    # pysam requires at least 2 alleles
    if len(new_alleles) < 2:   
        assert len(new_alleles) == 1
        for a in "ACGTN":
            if a != new_alleles[0]:
                break
        new_alleles.append(a)

    # new_record() refers to 
    # https://github.com/pysam-developers/pysam/blob/c6275315a8086858dae67f94c5caaf79a020fc46/pysam/libcbcf.pyx#L2079
    new_rec = out.new_record(
        contig = rec.chrom,
        start = rec.start,
        stop = rec.stop,
        alleles = new_alleles,
        id = None,
        qual = None,
        filter = "PASS",
        info = {"OAL":",".join(rec.alleles)},
        samples = new_samples
    )
    for i in range(len(rec.samples)):
        new_rec.samples[i].phased = rec.samples[i].phased

    return((ret, new_rec))


def __fix_file(in_fn, out_fn, ref_fn, verbose = False):
    """Fix REF, ALT and GT while delete other fields in FORMAT

    Parameters
    ----------
    in_fn : str
        Input VCF file to be fixed.
    out_fn : str
        Output VCF file.
    ref_fn : str
        Reference genome fasta file.
    verbose : bool
        Whether to print detailed log info.

    Returns
    -------
    int
        0 if success, -1 otherwise.
    """
    # load ref fasta 
    FASTA = pysam.FastaFile(ref_fn)   # if .fai file doesnot exist, this command would create one automatically.

    # load input vcf 
    iv = pysam.VariantFile(in_fn, "r")

    # add headers to output vcf
    header = iv.header.copy()
    header.add_line('##INFO=<ID=OAL,Number=.,Type=String,Description="Old alleles in the order of REF,ALTs">')
    header.add_line('##FORMAT=<ID=OGT,Number=1,Type=String,Description="Old Genotype">')

    if not out_fn:
        out_fn = "-"      # stdout
    ov = pysam.VariantFile(out_fn, "w", header = header)
    
    # fixref each vcf record and output
    nr = 0
    valid_cnt = 0
    fix_cnt = 0
    matched_cnt = 0
    for rec in iv.fetch():
        nr += 1
        try:
            ref_allele = FASTA.fetch(rec.chrom, rec.pos - 1, rec.pos)
            if not ref_allele or ref_allele.upper() not in ('A', 'C', 'G', 'T', 'N'):
                raise ValueError
            ref_allele = ref_allele.upper()
        except (IndexError, ValueError):
            if verbose:
                sys.stderr.write("[W::%s] %s:%d-%d invalid contig or coordinates out of range\n" % (
                                 COMMAND, rec.chrom, rec.pos, rec.pos))
            continue
            
        ret, new_rec = __fix_rec(ref_allele, rec, ov, verbose)
        if ret < 0:
            continue
        else:
            valid_cnt += 1
            if ret == 0: matched_cnt += 1
            else: fix_cnt += 1
            ov.write(new_rec)
        
    iv.close()
    ov.close()

    sys.stderr.write("%d total records in input VCF!\n" % (nr))
    sys.stderr.write("%d valid records in input VCF!\n" % (valid_cnt))
    sys.stderr.write("%d records have been fixed REF!\n" % (fix_cnt))
    sys.stderr.write("%d records don't need to fix REF!\n" % (matched_cnt))

    return(0)


def __usage(fp = sys.stdout):
    msg =  "\n"
    msg += "Version: %s\n" % VERSION
    msg += "Usage:   %s %s [options]\n" % (APP, COMMAND)
    msg += "\n"                                                        \
           "Options:\n"                                                \
           "  -i, --input FILE    Path to input vcf file\n"              \
           "  -r, --ref FILE      Path to reference genome fasta file\n"     \
           "  -o, --output FILE   Path to output vcf file. if not set, output to stdout\n"    \
           "  -v, --verbose       If use, output more detailed log info\n"                   \
           "  -h, --help          Print this message\n"                                       \
           "\n"                                            \
           "Note:\n"                                       \
           "1. REF, ALT and GT would be checked and possibly fixed. Fields other than GT\n" \
           "   in FORMAT would be deleted.\n"                     \
           "\n"    
    fp.write(msg)


def fixref_main(argv):
    # parse and check command line
    if len(argv) < 3:
        __usage(sys.stdout)
        sys.exit(0)
           
    opts, args = getopt.getopt(
        argv[2:], 
        "-h-i:-r:-o:-v", 
        ["help", "input=", "ref=", "output=", "verbose"])
    ref_file = in_vcf_file = out_vcf_file = None
    verbose = False
    for op, val in opts:
        if op in ("-i", "--input"): in_vcf_file = val
        elif op in ("-r", "--ref"): ref_file = val
        elif op in ("-o", "--output"): out_vcf_file = val
        elif op in ("-v", "--verbose"): verbose = True
        elif op in ("-h", "--help"): __usage(sys.stdout); sys.exit(0)
        else: sys.stderr.write("invalid option: %s\n" % op); sys.exit(1)

    # TODO: check args
    
    __fix_file(in_vcf_file, out_vcf_file, ref_file, verbose)