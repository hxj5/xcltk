# utils.py - help functions.

import os
import sys

from ...utils.grange import Region, RegionSet
from ...utils.zfile import zopen, ZF_F_GZIP, ZF_F_PLAIN


def load_region_from_txt(fn, sep = "\t", verbose = False):
    """Load regions from plain file.

    Parameters
    ----------
    fn : str
        Path to header-free plain file listing regions, each per line.
        The first 4 columns should be
        <chrom>, <start>, <end> (both 1-based, inclusive), <name>.
    verbose : bool
        Whether to print detailed log info.

    Returns
    -------
    list
        A list of `BlockRegion` objects if success, `None` otherwise.    
    """
    func = "load_region_from_txt"
    fp = zopen(fn, "rt")
    reg_list = []
    nl = 0
    if verbose:
        sys.stderr.write("[I::%s] start to load regions from file '%s' ...\n" % (func, fn))
    for line in fp:
        nl += 1
        parts = line.rstrip().split(sep)
        if len(parts) < 4:
            if verbose:
                sys.stderr.write("[E::%s] too few columns of line %d.\n" % (func, nl))
            return None           
        chrom, start, end, name = parts[:4]
        start, end = int(start), int(end)
        reg = Region(chrom, start, end + 1, name)
        reg_list.append(reg)
    fp.close()
    return reg_list


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