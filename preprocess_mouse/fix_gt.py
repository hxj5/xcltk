#!/usr/bin/env python

import sys
import getopt
import time
import gzip
import pysam


def __fix_rec(rec, gt0, phased, out, verbose = False):
    """
    @abstract       Check & fix GT for one vcf record
    @param rec      The vcf record whose ref is to be checked [pysam.VariantRecord]
    @param gt0      Genotype [tuple<int, int>]
    @param phased   Whether the GTs have been phased [bool]
    @param out      The output vcf object [pysam.VariantFile]
    @param verbose  Whether to print detailed log info [bool]
    @return         A tuple of two elements [tuple<int, pysam.VariantRecord>]: 
                      the running state and 
                      the checked vcf record 
                    The running state: 
                      -1, if error;
                      0, if the vcf record is successfully fixed.
    """
    new_samples = []
    gt = None
    for sample in rec.samples.values():
        if "GT" in rec.format: 
            sep = "|" if sample.phased else "/"
            ra_idx = sample.get("GT", ())      # index of ref/alt, a tuple
            gt = sep.join([str(a) for a in ra_idx if a is not None])
            if not gt: gt = "."
        else:
            gt = "."
        new_samples.append({
            "OGT":gt,
            "GT":gt0
        })

    # new_record() refers to 
    # https://github.com/pysam-developers/pysam/blob/c6275315a8086858dae67f94c5caaf79a020fc46/pysam/libcbcf.pyx#L2079
    new_rec = out.new_record(
        contig = rec.chrom,
        start = rec.start,
        stop = rec.stop,
        alleles = rec.alleles,
        id = None,
        qual = None,
        filter = "PASS",
        info = {},
        samples = new_samples
    )
    for i in range(len(rec.samples)):
        new_rec.samples[i].phased = phased

    ret = 0
    return((ret, new_rec))


def __fix_file(in_fn, out_fn, gt0, phased = False, verbose = False):
    """
    @abstract       Fix GT.
    @param in_fn    Input vcf file to be fixed [str]
    @param out_fn   Output vcf file [str]
    @param gt0      Genotype [tuple<int, int>]
    @param phased   Whether the GTs have been phased [bool]
    @param verbose  Whether to print detailed log info [bool]
    @return         0 if success, -1 otherwise
    """
    # load input vcf 
    iv = pysam.VariantFile(in_fn, "r")

    # add headers to output vcf
    header = iv.header.copy()
    header.add_line('##FORMAT=<ID=OGT,Number=1,Type=String,Description="Old Genotype">')

    if not out_fn:
        out_fn = "-"      # stdout
    ov = pysam.VariantFile(out_fn, "w", header = header)
    
    # fixref each vcf record and output
    nr = 0
    valid_cnt = 0
    for rec in iv.fetch():
        nr += 1
        ret, new_rec = __fix_rec(rec, gt0, phased, ov, verbose)
        if ret < 0:
            continue
        else:
            valid_cnt += 1
            ov.write(new_rec)
        
    iv.close()
    ov.close()

    sys.stderr.write("%d total records in input VCF!\n" % (nr, ))
    sys.stderr.write("%d valid records in input VCF!\n" % (valid_cnt, ))

    return(0)


def __usage(fp = sys.stderr):
    msg =  "\n"
    msg += "Usage: %s [options]\n" % (APP, )
    msg += "\n"                                                        \
           "Options:\n"                                                \
           "  -i, --input FILE    Path to input vcf file\n"              \
           "  -g, --gt STR        Genotype\n"     \
           "  -o, --output FILE   Path to output vcf file. if not set, output to stdout\n"    \
           "  -v, --verbose       If use, output more detailed log info\n"                   \
           "  -h, --help          Print this message\n"                                       \
           "\n"
    fp.write(msg)


def fix_gt(argv):
    # parse and check command line
    if len(argv) < 3:
        __usage(sys.stderr)
        sys.exit(1)
           
    opts, args = getopt.getopt(argv[1:], "-h-i:-g:-o:-v", ["help", "input=", "gt=", "output=", "verbose"])
    gt = in_vcf_file = out_vcf_file = None
    verbose = False
    for op, val in opts:
        if op in ("-i", "--input"): in_vcf_file = val
        elif op in ("-g", "--gt"): gt = val
        elif op in ("-o", "--output"): out_vcf_file = val
        elif op in ("-v", "--verbose"): verbose = True
        elif op in ("-h", "--help"): __usage(sys.stderr); sys.exit(1)
        else: sys.stderr.write("invalid option: %s\n" % op); sys.exit(1)

    # TODO: check args
    
    sep = "/" if "/" in gt else "|"
    gt0 = tuple([int(i) for i in gt.split(sep)])
    phased = sep == "|"
    __fix_file(in_vcf_file, out_vcf_file, gt0, phased, verbose)


APP = "fix_gt.py"


if __name__ == "__main__":
    fix_gt(sys.argv)

