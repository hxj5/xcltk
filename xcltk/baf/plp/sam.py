# Sam File Routine
# Author: Xianjie Huang

def check_read(read, conf):
    if read.mapq < conf.min_mapq:
        return(-2)
    if conf.excl_flag and read.flag & conf.excl_flag:
        return(-3)
    if conf.incl_flag and not read.flag & conf.incl_flag:
        return(-4)
    if conf.no_orphan and read.flag & BAM_FPAIRED and not \
        read.flag & BAM_FPROPER_PAIR:
        return(-5)
    if conf.cell_tag and not read.has_tag(conf.cell_tag):
        return(-11)
    if conf.umi_tag and not read.has_tag(conf.umi_tag):
        return(-12)
    if len(read.positions) < conf.min_len:
        return(-21)
    return(0)

def sam_fetch(sam, chrom, start, end):
    """Provide a wrapper for sam-fetch method that could automatically
       handle chrom with or without "chr" prefix.
    @param sam    A pysam.AlignmentFile object.
    @param chrom  Chromosome name [str]
    @param start  1-based, inclusive [int]
    @param end    1-based, inclusive [int]
    @return       Iterator if success, None otherwise. 
    """
    try:   # sam.fetch(): start and stop denote 0-based, half-open intervals.
        itr = sam.fetch(chrom, start - 1, end) 
    except:
        pass
    else:
        if itr:
            return itr
    chrom = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
    try:
        itr = sam.fetch(chrom, start - 1, end)
    except:
        return None
    else:
        return itr if itr else None

BAM_FPAIRED = 1
BAM_FPROPER_PAIR = 2

# Cigar
# reference: https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples
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
