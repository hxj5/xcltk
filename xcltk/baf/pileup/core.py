# core.py

import math
import os
import pickle
import pysam
import sys

from .mcount import MCount, MCNT_F_SPLICED, MCNT_F_UNSPLICED, MCNT_F_AMBIGUOUS
from .sam import check_read, sam_fetch
from .utils import load_jcset_from_gff, load_jcset_from_sam
from .zfile import zopen, ZF_F_GZIP

def sp_gene(thdata, mcnt, conf, gene_idx, fp_jc, fp_reg, fp_sp, fp_us, fp_am):
    gene = mcnt.gene
    for sam in conf.sam_list:
        itr = sam_fetch(sam, gene.chrom, gene.start, gene.end)
        if not itr:     # TODO: still need to output gene info to file.
            return(-3)
        for read in itr:
            if check_read(read, conf) < 0:
                continue
            ret = mcnt.push_read(read)
            if ret < 0:
                if ret == -1:
                    return(-5)
                continue

    if mcnt.stat() < 0:
        return(-7)

    if mcnt.jdata is not None and mcnt.jdata.size > 0:
        s = ""
        jc_list = mcnt.intv.get_junctions(sort = True)
        for jc in jc_list:
            n_jumi, n_jcell, is_valid = mcnt.jdata[jc.index, :]
            s += "\t".join([jc.chrom, str(jc.start), str(jc.end), 
                            jc.tran_id, jc.gene_id, 
                            str(n_jumi), str(n_jcell), str(is_valid)]) + "\n"
        fp_jc.write(s)

    gene = mcnt.gene
    nu_sp, nu_us, nu_am = [mcnt.tcount[x] for x in 
        (MCNT_F_SPLICED, MCNT_F_UNSPLICED, MCNT_F_AMBIGUOUS)]
    s = "\t".join([gene.chrom, str(gene.start), str(gene.end), 
                   gene.gene_id, gene.gene_name, 
                   str(nu_sp), str(nu_us), str(nu_am)]) + "\n"
    fp_reg.write(s)

    str_sp, str_us, str_am = "", "", ""
    for i, smp in enumerate(conf.barcodes):
        scnt = mcnt.cell_cnt[smp]
        nu_sp, nu_us, nu_am = [scnt.tcount[x] for x in 
            (MCNT_F_SPLICED, MCNT_F_UNSPLICED, MCNT_F_AMBIGUOUS)]
        if nu_sp > 0:
            str_sp += "%d\t%d\t%d\n" % (gene_idx, i + 1, nu_sp)
            thdata.nr_sp += 1
        if nu_us > 0:
            str_us += "%d\t%d\t%d\n" % (gene_idx, i + 1, nu_us)
            thdata.nr_us += 1
        if nu_am > 0:
            str_am += "%d\t%d\t%d\n" % (gene_idx, i + 1, nu_am)
            thdata.nr_am += 1
    fp_sp.write(str_sp)
    fp_us.write(str_us)
    fp_am.write(str_am)
    
    return(0)

# TODO: use clever IPC (Inter-process communication) instead of naive `raise Error`.
# NOTE: 
# 1. in current version, the gene idx is based on all genes (no gene filtering 
#    such as --minCOUNT)
# 2. bgzf errors when using pysam.AlignmentFile.fetch in parallel (with multiprocessing)
#    https://github.com/pysam-developers/pysam/issues/397
def sp_count(thdata):
    func = "sp_count"
    conf = thdata.conf
    thdata.ret = -1

    conf.sam_list = []
    for sam_fn in conf.sam_fn_list:
        sam = pysam.AlignmentFile(sam_fn, "r")    # auto detect file format
        conf.sam_list.append(sam)

    reg_list = None
    if thdata.is_reg_pickle:
        with open(thdata.reg_obj, "rb") as fp:
            reg_list = pickle.load(fp)
        os.remove(thdata.reg_obj)
    else:
        reg_list = thdata.reg_obj
    thdata.nr_reg = len(reg_list)

    fp_jc = zopen(thdata.out_junction_fn, "wt", ZF_F_GZIP, is_bytes = False)
    fp_reg = zopen(thdata.out_region_fn, "wt", ZF_F_GZIP, is_bytes = False)
    fp_sp = zopen(thdata.out_spliced_fn, "wt", ZF_F_GZIP, is_bytes = False)
    fp_us = zopen(thdata.out_unspliced_fn, "wt", ZF_F_GZIP, is_bytes = False)
    fp_am = zopen(thdata.out_ambiguous_fn, "wt", ZF_F_GZIP, is_bytes = False)

    mcnt = MCount(conf.barcodes, conf)
    m_reg = float(len(reg_list))
    n_reg = 0
    l_reg = 0
    min_jc = None
    for gene_idx, gene in enumerate(reg_list):
        if conf.debug > 0:
            sys.stderr.write("[D::%s][Thread-%d] processing gene '%s' ...\n" %
                              (func, thdata.idx, gene.gene_id))

        rs = load_jcset_from_gff(gene)
        if rs is None:
            raise ValueError("[%s] errcode %d" % (func, -3))
        if conf.debug > 0:
            sys.stderr.write("[D::%s][Thread-%d] gene %s: %d junctions extracted\n" %
                              (func, thdata.idx, gene.gene_id, rs.get_n()))
            if conf.debug > 1:
                junctions = rs.get_junctions(sort = True)
                for jc in junctions:
                    sys.stderr.write("\t%s\t%s\t%s\t%s\t%s\n" % 
                        (jc.index, jc.chrom, jc.start, jc.end, jc.tran_id))

        if conf.incl_denovo:  # CHECK ME! how STAR detect and represent junctions.
            if conf.denovo_jc_len >= 1:
                min_jc = conf.denovo_jc_len
            elif conf.denovo_jc_len > 0:
                min_jc = rs.get_len_quantile(conf.denovo_jc_len)
            else:    # <= 0 means no filtering for length of denovo junctions.
                min_jc = conf.denovo_jc_len
            for sam in conf.sam_list:
                rs1 = load_jcset_from_sam(sam, gene.chrom, gene.start, gene.end, conf, min_jc)
                if rs1 is None:
                    raise ValueError("[%s] errcode %d" % (func, -5))
                rs.merge(rs1)
            if conf.debug > 0:
                sys.stderr.write("[D::%s][Thread-%d] gene %s: %d junctions including denovo ones\n" %
                                  (func, thdata.idx, gene.gene_id, rs.get_n()))

        if mcnt.add_gene(gene, rs) < 0:
            raise ValueError("[%s] errcode %d" % (func, -7))
            
        if sp_gene(thdata, mcnt, conf, gene_idx + 1, fp_jc, fp_reg, fp_sp, fp_us, fp_am) < 0:
            raise ValueError("[%s] errcode %d" % (func, -9))
        mcnt.reset()

        n_reg += 1
        frac_reg = n_reg / m_reg
        if frac_reg - l_reg >= 0.02 or n_reg == m_reg:
            sys.stdout.write("[I::%s][Thread-%d] %d%% genes processed\n" % 
                (func, thdata.idx, math.floor(frac_reg * 100)))
            l_reg = frac_reg

    fp_jc.close()
    fp_reg.close()
    fp_sp.close()
    fp_us.close()
    fp_am.close()

    thdata.conf = None    # sam object cannot be pickled.
    thdata.ret = 0

    if thdata.out_fn:
        with open(thdata.out_fn, "wb") as fp_td:
            pickle.dump(thdata, fp_td)
            
    return((0, thdata))
