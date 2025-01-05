# core.py - core part of feature counting.

import math
import os
import pickle
import pysam

from logging import debug, info

from .mcount import MCount
from ...utils.sam import sam_fetch, \
    BAM_FPAIRED, BAM_FPROPER_PAIR
from ...utils.zfile import zopen, ZF_F_GZIP


def __get_include_len_given_range(r1, r2):
    # get the length of included part of r1 within r2.
    # we already know that r1 overlaps with r2.
    s1, e1 = r1[:2]
    s2, e2 = r2[:2]
    if s1 < s2 and e1 > e2:
        return(0)
    if s1 < s2:
        assert e1 >= s2
        return(e1 - s2)
    if e1 > e2:
        assert s1 <= e2
        return(e2 - s1)
    return(e1 - s1)


def __get_include_frac(pos_list, s, e):
    n = len(pos_list)
    if n <= 0:
        return(None)
    m = __get_include_len(pos_list, s, e)
    return(m / float(n))


def __get_include_len(pos_list, s, e):
    # all input parameters are 0-based.
    include_pos_list = [x for x in pos_list if s <= x <= e]
    return(len(include_pos_list))


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


# TODO: use clever IPC (Inter-process communication) instead of naive `raise Error`.
# NOTE: 
# 1. bgzf errors when using pysam.AlignmentFile.fetch in parallel (with multiprocessing)
#    https://github.com/pysam-developers/pysam/issues/397
def fc_features(thdata):
    conf = thdata.conf
    thdata.ret = -1

    sam_list = []
    for sam_fn in conf.sam_fn_list:
        sam = pysam.AlignmentFile(sam_fn, "r")    # auto detect file format
        sam_list.append(sam)

    reg_list = None
    if thdata.is_reg_pickle:
        with open(thdata.reg_obj, "rb") as fp:
            reg_list = pickle.load(fp)
        os.remove(thdata.reg_obj)
    else:
        reg_list = thdata.reg_obj

    fp_reg = zopen(thdata.out_region_fn, "wt", ZF_F_GZIP, is_bytes = False)
    fp_mtx = zopen(thdata.out_mtx_fn, "wt", ZF_F_GZIP, is_bytes = False)

    mcnt = MCount(conf.samples, conf)

    m_reg = float(len(reg_list))
    n_reg = 0         # number of processed genes.
    l_reg = 0         # fraction of processed genes, used for verbose.
    k_reg = 1         # index of output region in sparse matrix, 1-based.

    for reg_idx, reg in enumerate(reg_list):
        if conf.debug > 0:
            debug("[Thread-%d] processing region '%s' ..." % \
                (thdata.idx, reg.get_id()))
            
        mcnt.reset()
        mcnt.add_region(reg)
        str_reg = "%s\t%d\t%d\t%s\n" % \
            (reg.chrom, reg.start, reg.end - 1, reg.get_id())
        ret, counts = fc_fet1(reg, sam_list, mcnt, conf)
        if ret < 0:
            raise RuntimeError("errcode -9")

        str_mtx = ""
        for i, smp in enumerate(conf.samples):
            nu = counts[smp]
            if nu <= 0:
                continue
            else:
                str_mtx += "%d\t%d\t%d\n" % (k_reg, i + 1, nu)
                thdata.nr_mtx += 1

        if str_mtx:
            fp_mtx.write(str_mtx)
            fp_reg.write(str_reg)
            k_reg += 1
        elif conf.output_all_reg:
            fp_reg.write(str_reg)
            k_reg += 1

        n_reg += 1
        frac_reg = n_reg / m_reg
        if frac_reg - l_reg >= 0.02 or n_reg == m_reg:
            info("[Thread-%d] %d%% genes processed" % 
                (thdata.idx, math.floor(frac_reg * 100)))
            l_reg = frac_reg

    thdata.nr_reg = k_reg - 1

    fp_reg.close()
    fp_mtx.close()
    for sam in sam_list:
        sam.close()
    sam_list.clear()

    thdata.conf = None    # sam object cannot be pickled.
    thdata.ret = 0

    if thdata.out_fn:
        with open(thdata.out_fn, "wb") as fp_td:
            pickle.dump(thdata, fp_td)
            
    return((0, thdata))


def fc_fet1(reg, sam_list, mcnt, conf):
    ret = None
    for idx, sam in enumerate(sam_list):
        itr = sam_fetch(sam, reg.chrom, reg.start, reg.end - 1)
        if not itr:    
            continue
        for read in itr:
            if check_read(read, conf) < 0:
                continue
            if 0 < conf.min_include < 1:
                if __get_include_frac(read.positions, reg.start - 1, reg.end - 2) < conf.min_include:
                    continue
            else:
                if __get_include_len(read.positions, reg.start - 1, reg.end - 2) < conf.min_include:
                    continue
            if conf.use_barcodes():
                ret = mcnt.push_read(read)
            else:
                sample = conf.samples[idx]
                ret = mcnt.push_read(read, sample)
            if ret < 0:
                if ret == -1:
                    return((-5, None))
                continue
    if mcnt.stat() < 0:
        return((-7, None))
    counts = {smp:mcnt.cell_cnt[smp].tcount for smp in conf.samples}
    return((0, counts))
