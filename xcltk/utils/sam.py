# Copied from cellSNP, https://github.com/single-cell-genetics/cellSNP/blob/purePython/cellSNP/utils/cellsnp_utils.py
# SAM file utils 
# Author: Xianjie Huang

import pysam

# Operations of CIGAR, copied from pysam (https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples)
# Do not change these values unless you know what you are doing.
BAM_CMATCH = 0
BAM_CINS = 1
BAM_CDEL = 2
BAM_CREF_SKIP = 3
BAM_CSOFT_CLIP = 4
BAM_CHARD_CLIP = 5
BAM_CPAD = 6
BAM_CEQUAL = 7
BAM_CDIFF = 8
BAM_CBACK = 9


# for pysam
CACHE_CHROM = None
CACHE_SAMFILE = None

def check_pysam_chrom(samFile, chrom=None):
    """Chech if samFile is a file name or pysam object, and if chrom format. 
    """
    global CACHE_CHROM
    global CACHE_SAMFILE

    if CACHE_CHROM is not None:
        if (samFile == CACHE_SAMFILE) and (chrom == CACHE_CHROM):
            return CACHE_SAMFILE, CACHE_CHROM

    if type(samFile) == str:
        ftype = samFile.split(".")[-1]
        if ftype != "bam" and ftype != "sam" and ftype != "cram" :
            print("Error: file type need suffix of bam, sam or cram.")
            sys.exit(1)
        if ftype == "cram":
            samFile = pysam.AlignmentFile(samFile, "rc")
        elif ftype == "bam":
            samFile = pysam.AlignmentFile(samFile, "rb")
        else:
            samFile = pysam.AlignmentFile(samFile, "r")

    if chrom is not None:
        if chrom not in samFile.references:
            if chrom.startswith("chr"):
                chrom = chrom.split("chr")[1]
            else:
                chrom = "chr" + chrom
        if chrom not in samFile.references:
            print("Can't find references %s in samFile" %chrom)
            return samFile, None
    
    CACHE_CHROM = chrom
    CACHE_SAMFILE = samFile
    return samFile, chrom

def get_query_bases(read, full_length=False):
    """
    @abstract            Return a list of bases in qurey sequence that are within the alignment.
    @param read          An AlignedSegment object. [AlignedSegment]
    @param full_length   If full_length is set, None values will be included for any soft-clipped or 
                         unaligned positions within the read. The returned list will thus be of the 
                         same length as the read. [bint]
    @return              A list of bases. [list]
    """
    cigar_tuples = read.cigartuples
    if not cigar_tuples:
        return []

    result = []
    pos = 0
    s = read.query_sequence

    for op, l in cigar_tuples:
        if op == BAM_CSOFT_CLIP or op == BAM_CINS:
            if full_length:
                for i in range(0, l):
                    result.append(None)
            pos += l
        elif op == BAM_CMATCH or op == BAM_CEQUAL or op == BAM_CDIFF:
            for i in range(pos, pos + l):
                result.append(s[i])
            pos += l
        # else: do nothing.
    return result

def get_query_qualities(read, full_length=False):
    """
    @abstract            Return a list of qualities in qurey quality sequence that are within the alignment.
    @param read          An AlignedSegment object. [AlignedSegment]
    @param full_length   If full_length is set, None values will be included for any soft-clipped or 
                         unaligned positions within the read. The returned list will thus be of the 
                         same length as the read. [bint]
    @return              A list of qualities. [list]

    @note                The returned qual values are not ASCII-encoded values typically seen in FASTQ or SAM formatted files,
                         so no need to substract 33.
    """
    cigar_tuples = read.cigartuples
    if not cigar_tuples:
        return []

    result = []
    pos = 0
    s = read.query_qualities

    for op, l in cigar_tuples:
        if op == BAM_CSOFT_CLIP or op == BAM_CINS:
            if full_length:
                for i in range(0, l):
                    result.append(None)
            pos += l
        elif op == BAM_CMATCH or op == BAM_CEQUAL or op == BAM_CDIFF:
            for i in range(pos, pos + l):
                result.append(s[i])
            pos += l
        # else: do nothing.
    return result

