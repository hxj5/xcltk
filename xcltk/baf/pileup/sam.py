# Sam File Routine
# Author: Xianjie Huang

class CigarBlock:
    """Cigar Block
    Parameters
    ----------
    op: cigar operation [int]     
    length: length of the block (i.e., the number in the cigar) [int]
    start: genomic coordinate, 1-based, exclusive [int]
    end: genomic coordinate, 1-based, exclusive [int]  
    qstart: query index, 0-based, exclusive [int]
    qend: query index, 0-based, exclusive [int]
    qseq: query read base sequence; "" for DEL/REF_SKIP [str]
    qqual: query read quality sequence; "" for DEL/REF_SKIP [str]
    """
    def __init__(self, op, length, start, end, 
                 qstart, qend, qseq, qqual):
        self.op = op
        self.length = length
        self.start = start  
        self.end = end 
        self.qstart = qstart
        self.qend = qend  
        self.qseq = qseq  
        self.qqual = qqual

def extract_cigar_blocks(read):
    """Parse one read
    @param read   A pysam::AlignedSegment object.
    @return       A list of CigarBlock objects.
    """
    cigar_tuples = read.cigartuples
    if not cigar_tuples:
        return []

    result = []
    pos = read.reference_start    # left most genomic coordinate, 0-based
    qpos = 0
    qseq = read.query_sequence
    qqual = read.query_qualities

    for op, l in cigar_tuples:
        # see class CigarBlock for the details of each c_xxx variable.
        c_qstart = c_qend = -1
        c_qseq = c_qqual = ""
        c_start = c_end = -1
        if op in (BAM_CSOFT_CLIP, BAM_CINS):
            c_qstart = qpos - 1
            c_qend = qpos + l
            c_qseq = qseq[qpos:(qpos + l)]
            c_qqual = qqual[qpos:(qpos + l)]
            c_start = pos
            c_end = pos + 1
            qpos += l
        elif op in (BAM_CMATCH, BAM_CEQUAL, BAM_CDIFF):
            c_qstart = qpos - 1
            c_qend = qpos + l
            c_qseq = qseq[qpos:(qpos + l)]
            c_qqual = qqual[qpos:(qpos + l)]
            c_start = pos
            c_end = pos + l + 1
            pos += l
            qpos += l
        elif op in (BAM_CDEL, BAM_CREF_SKIP):
            c_qstart = qpos - 1
            c_qend = qpos
            c_start = pos
            c_end = pos + l + 1
            pos += l
        #else: do nothing
        cigar_block = CigarBlock(
            op, l,
            c_start, c_end,
            c_qstart, c_qend, c_qseq, c_qqual
        )
        result.append(cigar_block)

    return result

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
