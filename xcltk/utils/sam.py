# sam.py - SAM/BAM file routine


def get_query_bases(read, full_length=False):
    """Qurey bases that are within the alignment.

    Parameters
    ----------
    read : pysam::AlignedSegment object
        The alignment to be queried.
    full_length : bool
        If full_length is set, `None` values will be included for any 
        soft-clipped or unaligned positions within the read. 
        The returned list will thus be of the same length as the read.
    
    Returns
    -------
    list
        A list of bases in qurey sequence that are within the alignment.
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
    """Qurey qualities that are within the alignment.

    Parameters
    ----------
    read : pysam::AlignedSegment object
        The alignment to be queried.
    full_length : bool
        If full_length is set, `None` values will be included for any 
        soft-clipped or unaligned positions within the read. 
        The returned list will thus be of the same length as the read.
    
    Returns
    -------
    list
        A list of bases in qurey sequence that are within the alignment.
        Note that the returned qual values are not ASCII-encoded values 
        typically seen in FASTQ or SAM formatted files, no need to 
        substract 33.
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


def sam_fetch(sam, chrom, start, end):
    """Provide a wrapper for sam-fetch method that could automatically
       handle chrom with or without "chr" prefix.

    Parameters
    ----------
    sam : pysam.AlignmentFile object.
        The BAM file object.
    chrom : str
        Chromosome name.
    start : int
        1-based, inclusive.
    end : int
        1-based, inclusive.

    Returns
    -------
    Iterator
        Iterator if success, `None` otherwise.
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
