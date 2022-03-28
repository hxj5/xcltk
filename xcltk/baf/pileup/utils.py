# Utils

import os
import sys
from .region import Junction, JunctionSet
from .sam import check_read, sam_fetch
from .zfile import zopen, ZF_F_GZIP, ZF_F_PLAIN

#TODO: implement the following two functions
#  note the format of the STAR output junction file.
#def load_jcset_from_bed()
#def load_jcset_from_tsv()

def load_jcset_from_gff(gene, verbose = True):   # @param verbose  Verbose level: 0, 1, 2 [int]
    """Load JunctionSet from a gff gene object
    @param gene     A gff::Gene object
    @param verbose  Verbose or not [bool]
    @return      A JunctionSet object if success, None otherwise.
    @note Duplicate junctions (e.g., from different transcripts) will
        be discarded.
    """
    func = "load_jcset_from_gff"
    rs = JunctionSet(is_uniq = True)
    chrom = gene.chrom
    for tran in gene.trans:
        for i in range(1, tran.exons.shape[0]):
            start = tran.exons[i - 1, 1] + 1
            end = tran.exons[i, 0]
            junc_id = "%s-%s-%s" % (chrom, start, end)
            if start >= end:
                if verbose:
                    sys.stderr.write("[W::%s] skip invalid junction %s.\n" % 
                        (func, junc_id))
                continue
            #<DEV/>#exon5 = (tran.exons[i - 1, 0], tran.exons[i - 1, 1] + 1)
            #<DEV/>#exon3 = (tran.exons[i, 0], tran.exons[i, 1] + 1)
            jc = Junction(
                chrom, start, end, 
                junc_id = junc_id, 
                tran_id = tran.tran_id, 
                gene_id = gene.gene_id,
                is_anno = True)
            rs.add(jc)
    return rs

def load_jcset_from_sam(sam, chrom, start, end, conf, min_jc):
    """Load junction set of specific region from sam file.
    @param sam     A pysam.AlignmentFile object.
    @param chrom   Chromosome name [str]
    @param start   1-based, inclusive [int]
    @param end     1-based, inclusive [int]
    @param conf    A config::Config object.
    @param min_jc  Minimum junction length [float]
    @return        A JunctionSet object if success, None otherwise.
    """
    func = "load_jcset_from_sam"
    itr = sam_fetch(sam, chrom, start, end)
    if not itr:
        return None
    
    rs = JunctionSet(is_uniq = True)
    for read in itr:
        if check_read(read, conf) < 0:
            continue
        cigar_blocks = extract_cigar_blocks(read)
        if not cigar_blocks:
            continue      
        for cb in cigar_blocks:
            if cb.op not in (BAM_CDEL, BAM_CREF_SKIP):
                continue
            if cb.length < min_jc:
                continue
            jstart = cb.start + 1
            jend = cb.end                  
            junc_id = "%s-%s-%s" % (chrom, jstart, jend)
            if jstart >= jend:
                if verbose:
                    sys.stderr.write("[W::%s] skip invalid junction %s.\n" % 
                        (func, junc_id))
                continue
            jc = Junction(
                chrom, jstart, jend, 
                junc_id = junc_id,
                tran_id = "", 
                gene_id = "",
                is_anno = False)
            rs.add(jc)
    return(rs)

def _fmt_line(ln, k):
        items = ln.split("\t")
        items[0] = str(int(items[0]) + k)
        return("\t".join(items))

# internal use only!
def merge_mtx(in_fn_list, in_format,
              out_fn, out_fmode, out_format,
              nrow_list, ncol, nrecord, remove = False):
    if len(in_fn_list) != len(nrow_list):
        return(-1)
    bufsize = 1048576    # 1M

    is_bytes = "b" in out_fmode
    out_fp = zopen(out_fn, out_fmode, out_format, is_bytes = is_bytes)

    nrow_total = sum(nrow_list)
    header  = "%%MatrixMarket matrix coordinate integer general\n"
    header += "%%\n"
    header += "%d\t%d\t%d\n" % (nrow_total, ncol, nrecord)
    if is_bytes:
        header = bytes(header, "utf8")
    out_fp.write(header)

    nline = 0
    k = 0
    for in_fn, nrow in zip(in_fn_list, nrow_list):
        with zopen(in_fn, "rt", in_format) as in_fp:
            while True:
                lines = in_fp.readlines(bufsize)
                if not lines:
                    break
                nline += len(lines)
                lines = [_fmt_line(ln, k) for ln in lines]
                s = "".join(lines)
                if is_bytes:
                    s = bytes(s, "utf8")
                out_fp.write(s)
        k += nrow
    out_fp.close()
    if nline != nrecord:
        return(-1)
    if remove:
        for in_fn in in_fn_list:
            os.remove(in_fn)
    return(0) 

# internal use only!
def merge_tsv(in_fn_list, in_format, 
              out_fn, out_fmode, out_format, 
              remove = False):
    bufsize = 1048576   # 1M
    is_bytes = "b" in out_fmode
    out_fp = zopen(out_fn, out_fmode, out_format, is_bytes)
    in_fmode = "rb" if is_bytes else "rt"
    for in_fn in in_fn_list:
        with zopen(in_fn, in_fmode, in_format) as in_fp:
            while True:
                dat = in_fp.read(bufsize)
                if not dat:
                    break
                out_fp.write(dat)
    out_fp.close()
    if remove:
        for in_fn in in_fn_list:
            os.remove(in_fn)
    return(0)

# internal use only!
def rewrite_mtx(in_fn, in_format, 
                out_fn, out_fmode, out_format, 
                nrow, ncol, nrecord, remove = False):
    if in_fn == out_fn:
        return(-1)
    bufsize = 1048576   # 1M
  
    is_bytes = "b" in out_fmode
    out_fp = zopen(out_fn, out_fmode, out_format, is_bytes = is_bytes)
    header  = "%%MatrixMarket matrix coordinate integer general\n"
    header += "%%\n"
    header += "%d\t%d\t%d\n" % (nrow, ncol, nrecord)
    if is_bytes:
        header = bytes(header, "utf8")
    out_fp.write(header)

    nline = 0
    in_fmode = "rb" if is_bytes else "rt"
    sep = b"" if is_bytes else ""
    with zopen(in_fn, in_fmode, in_format) as in_fp:
        while True:
            lines = in_fp.readlines(bufsize)
            if not lines:
                break
            nline += len(lines)
            s = sep.join(lines)
            out_fp.write(s)
    out_fp.close()
    if nline != nrecord:
        return(-1)
    if remove:
        os.remove(in_fn)
    return(0)
