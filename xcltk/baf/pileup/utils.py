# Utils

import os
import sys
from .region import Junction, JunctionSet
from .sam import check_read, sam_fetch
from .zfile import zopen, ZF_F_GZIP, ZF_F_PLAIN

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
