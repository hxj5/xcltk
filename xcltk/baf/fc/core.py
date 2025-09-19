# core.py - core part of feature counting.


import math
import os
import pickle
import pysam

from logging import debug, error, info

from .mcount import MCount
from ...utils.sam import sam_fetch, \
    BAM_FPAIRED, BAM_FPROPER_PAIR
from ...utils.zfile import zopen, ZF_F_GZIP



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
    fp_ad = zopen(thdata.out_ad_fn, "wt", ZF_F_GZIP, is_bytes = False)
    fp_dp = zopen(thdata.out_dp_fn, "wt", ZF_F_GZIP, is_bytes = False)
    fp_oth = zopen(thdata.out_oth_fn, "wt", ZF_F_GZIP, is_bytes = False)

    mcnt = MCount(conf.samples, conf)

    m_reg = float(len(reg_list))
    n_reg = 0         # number of processed genes.
    l_reg = 0         # fraction of processed genes, used for verbose.
    k_reg = 1         # index of output region in sparse matrix, 1-based.
    for reg_idx, reg in enumerate(reg_list):
        if conf.debug > 0:
            debug("[Thread-%d] processing region '%s' ..." % \
                (thdata.idx, reg.name))
            
        mcnt.reset()
        str_reg = "%s\t%d\t%d\t%s\n" % \
            (reg.chrom, reg.start, reg.end - 1, reg.name)
        if reg.snp_list:
            ret, reg_ref_cnt, reg_alt_cnt, reg_oth_cnt, reg_dp_cnt = \
                fc_fet1(reg, sam_list, mcnt, conf)
            if ret < 0:
                raise RuntimeError("errcode -9")

            str_ad, str_dp, str_oth = "", "", ""
            for i, smp in enumerate(conf.samples):
                nu_ad, nu_dp = reg_alt_cnt[smp], reg_dp_cnt[smp]
                nu_oth = reg_oth_cnt[smp]

                if nu_dp + nu_oth <= 0:
                    continue
                if nu_ad > 0:
                    str_ad += "%d\t%d\t%d\n" % (k_reg, i + 1, nu_ad)
                    thdata.nr_ad += 1
                if nu_dp > 0:
                    str_dp += "%d\t%d\t%d\n" % (k_reg, i + 1, nu_dp)
                    thdata.nr_dp += 1
                if nu_oth > 0:
                    str_oth += "%d\t%d\t%d\n" % (k_reg, i + 1, nu_oth)
                    thdata.nr_oth += 1

            if str_dp or str_oth:
                fp_ad.write(str_ad)
                fp_dp.write(str_dp)
                fp_oth.write(str_oth)
                fp_reg.write(str_reg)
                k_reg += 1
            elif conf.output_all_reg:
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
    fp_ad.close()
    fp_dp.close()
    fp_oth.close()
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
    reg_ref_umi = {smp:set() for smp in conf.samples}
    reg_alt_umi = {smp:set() for smp in conf.samples}
    reg_oth_umi = {smp:set() for smp in conf.samples}

    for snp in reg.snp_list:
        ret, mcnt = plp_snp(snp, sam_list, mcnt, conf)
        if ret < 0:
            error("SNP (%s:%d:%s:%s) pileup failed; errcode %d." % \
                (snp.chrom, snp.pos, snp.ref, snp.alt, ret))
            return((-3, None, None, None, None))
        elif ret > 0:     # snp filtered.
            continue
        for smp, scnt in mcnt.cell_cnt.items():
            for umi, ucnt in scnt.umi_cnt.items():
                if not ucnt.allele:
                    continue
                ale_idx = snp.get_region_allele_index(ucnt.allele)
                if ale_idx == 0:        # ref allele of the region.
                    reg_ref_umi[smp].add(umi)
                elif ale_idx == 1:      # alt allele of the region.
                    reg_alt_umi[smp].add(umi)
                else:
                    reg_oth_umi[smp].add(umi)

    reg_ref_cnt = {smp:0 for smp in conf.samples}
    reg_alt_cnt = {smp:0 for smp in conf.samples}
    reg_oth_cnt = {smp:0 for smp in conf.samples}
    reg_dp_cnt =  {smp:0 for smp in conf.samples}

    for smp in conf.samples:
        reg_ref_cnt[smp] = len(reg_ref_umi[smp])
        reg_alt_cnt[smp] = len(reg_alt_umi[smp])
        dp_umi = reg_ref_umi[smp].union(reg_alt_umi[smp]) # CHECK ME! theoretically no shared UMIs
        reg_oth_umi[smp] = reg_oth_umi[smp].difference(dp_umi)
        reg_oth_cnt[smp] = len(reg_oth_umi[smp])
        reg_dp_cnt[smp]  = len(dp_umi)

        if reg_ref_cnt[smp] + reg_alt_cnt[smp] != reg_dp_cnt[smp]:
            if conf.debug > 0:
                msg = "duplicate UMIs in region '%s', sample '%s':\n" % (
                    reg.name, smp)
                msg += "\tREF, ALT, DP_uniq (%d, %d, %d)." % (
                    reg_ref_cnt[smp], reg_alt_cnt[smp], reg_dp_cnt[smp])
                debug(msg)
            if conf.no_dup_hap:
                nu_share = reg_ref_cnt[smp] + reg_alt_cnt[smp] - reg_dp_cnt[smp]
                reg_ref_cnt[smp] = reg_ref_cnt[smp] - nu_share
                reg_alt_cnt[smp] = reg_alt_cnt[smp] - nu_share
            reg_dp_cnt[smp] = reg_ref_cnt[smp] + reg_alt_cnt[smp]
    
    return((0, reg_ref_cnt, reg_alt_cnt, reg_oth_cnt, reg_dp_cnt))



def plp_snp(snp, sam_list, mcnt, conf):
    """Pileup one SNP
    
    Parameters
    ----------
    snp : gfeature::SNP object
        The SNP to be pileuped.
    mcnt : mcount::MCount object
        The counting machine.
    conf : config::Config object
        Configuration.
    
    Returns
    -------
    tuple
        A tuple of 2 elements:
        (1) ret : int
            return code. 0 if success; negative if error; positive if filtered.
        (2) mcnt : mcount::MCount object
    """
    ret = None
    if mcnt.add_snp(snp) < 0:   # mcnt reset() inside.
        return((-3, mcnt))
    for idx, sam in enumerate(sam_list):
        itr = sam_fetch(sam, snp.chrom, snp.pos, snp.pos)
        if not itr:    
            continue
        for read in itr:
            if check_read(read, conf) < 0:
                continue
            if conf.use_barcodes():
                ret = mcnt.push_read(read)
            else:
                sample = conf.samples[idx]
                ret = mcnt.push_read(read, sample)
            if ret < 0:
                if ret == -1:
                    return((-5, mcnt))
                continue
    if mcnt.stat() < 0:
        return((-7, mcnt))
    snp_cnt = sum(mcnt.tcount)
    if snp_cnt < conf.min_count:
        return((3, mcnt))
    snp_ref_cnt = mcnt.tcount[mcnt.base_idx[snp.ref]]
    snp_alt_cnt = mcnt.tcount[mcnt.base_idx[snp.alt]]
    snp_minor_cnt = min(snp_ref_cnt, snp_alt_cnt)
    if snp_minor_cnt < snp_cnt * conf.min_maf:
        return((5, mcnt))
    return((0, mcnt))
