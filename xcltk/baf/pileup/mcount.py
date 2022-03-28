# Reads Counting
# Author: Xianjie Huang

# Note that the index of junctions should be double checked, especially
# in np.array (np.transpose)!

import sys
import numpy as np
from .mempool import MemPool
from .region import JunctionSet, format_chrom
from .sam import extract_cigar_blocks, check_read

class JCAlign:
    def __init__(self):  # this function takes no parameters as it's used by MemPool.
        self.sdata = None
        self.length = None
        self.gap = None

    def align(self, jc, read, positions, conf, overlap):
        """Align a set of genomic positions to one junction.
        @param jc        A region::Junction object.
        @param read      A pysam::AlignedSegment object.
        @param positions The aligned ref positions of @p read; 1-based [np.array]
        @param conf      A config::Config object.
        @param overlap   Does the positions overlap with @p jc [bool]
        @return          Void.
        """
        self.sdata = np.zeros((JC_R_N, ), dtype = ST_R_DTYPE)
       
        p5 = positions[positions < jc.start]
        pb = positions[(positions >= jc.start) & (positions < jc.end)]
        p3 = positions[positions >= jc.end]
        n5 = np.sum(p5 >= jc.exon5[0]) if jc.exon5 else p5.size  # CHECK ME! exon5/3 from various transcripts?
        nb = pb.size
        n3 = np.sum(p3 < jc.exon3[1]) if jc.exon3 else p3.size
        self.length = [n5, nb, n3]

        g5 = jc.start - 1 - np.max(p5) if p5.size > 0 else 0
        gb5 = np.min(pb) - jc.start if pb.size > 0 else 0
        gb3 = jc.end - 1 - np.max(pb) if pb.size > 0 else 0
        g3 = np.min(p3) - jc.end if p3.size > 0 else 0        
        self.gap = [g5, gb5, gb3, g3]

        if not overlap:
            return    # keep elements in stat all 0

        is_n5, is_nb, is_n3 = [int(x >= conf.min_overlap) for x in (n5, nb, n3)]
        is_g5, is_gb5, is_gb3, is_g3 = [int(x <= conf.max_gap) for x in (g5, gb5, gb3, g3)]

        if is_n5 + is_nb + is_n3 == 0 or is_n5 * (1 - is_nb) * (1 - is_n3) \
            or (1 - is_n5) * (1 - is_nb) * is_n3:
            return     # no overlapping
        self.sdata[JC_R_OVP] = 1

        is_p5 = is_n5 * is_g5
        is_pb5 = is_nb * is_gb5
        is_pb3 = is_nb * is_gb3
        is_p3 = is_n3 * is_g3

        if is_p5 * (1 - is_nb) * is_p3:
            self.sdata[JC_R_CROSS] = 1
            return
        if (1 - is_n5) * is_nb * (1 - is_n3):
            self.sdata[JC_R_IN] = 1
            return
        if is_p5 and is_pb5 :
            self.sdata[JC_R_BOUND5] = 1
        if is_pb3 and is_p3:
            self.sdata[JC_R_BOUND3] = 1

        if self.sdata[JC_R_BOUND5] + self.sdata[JC_R_BOUND3] == 0:
            self.sdata[JC_R_OTHER] = 1

    def reset(self):
        self.sdata = None
        self.length = None
        self.gap = None  

class JCReadAlign:
    def __init__(self):
        self.data = None

    def align(self, read, mcnt):
        """Align the read to a collection of junctions.
        @param read   A pysam::AlignedSegment object.
        @param mcnt   A MCount object containing the set of junctions.
        @return       0 if success, -1 otherwise [int]
        @note For now, we implement a naive alignment strategy, which
           1. does not incorportate the block type information (e.g., exon or intron).
           2. ...
        """
        func = "JCReadAlign::align"
        conf = mcnt.conf
        self.data = []  # a list of <jc, JCAlign> tuples.

        if conf.debug > 1:
            sys.stderr.write("[D::%s] aligning read '%s'.\n" %
                              (func, read.query_name))

        positions = np.array([x + 1 for x in read.get_reference_positions()])  # to 1-based.
        if len(positions) <= 0:
            return(0)

        junctions = mcnt.intv.fetch(
            read.reference_name,
            positions[0], 
            positions[-1] + 1)

        for jc in junctions:
            aln = mcnt.pool_aln.get()
            aln.align(jc, read, positions, mcnt.conf, True)
            if aln.sdata[JC_R_OVP]:
                self.data.append((jc, aln))

            if mcnt.conf.debug > 1:
                sys.stderr.write("\t[%s] sdata=%s; len=%s; gap=%s\n" % 
                                 (jc.get_id(), aln.sdata, aln.length, aln.gap))

        if mcnt.conf.debug > 1:
            sys.stderr.write("[D::%s] %d junctions overlaps with read.\n" % 
                              (func, len(self.data)))

        return(0)

    def reset(self):
        self.data = None

JC_R_OVP = 0         # whether overlap or not?
JC_R_CROSS = 1
JC_R_BOUND5 = 2
JC_R_IN = 3
JC_R_BOUND3 = 4
JC_R_OTHER = 5       # OVP and not any of [CROSS, B5, IN, B3]
JC_R_N = 6

class UCount:
    """Counting Unit
    @param scnt   A SCount object that the UCount object belongs to.
    @param conf   A config::Config object.
    @param intv   A region::JunctionSet object.
    @param pe     HashMap of <str:list> for reads_name:JCReadAlign pair for pair-end data.
    @param se     List of JCReadAlign objects for single-end data.
    @param n      Num of read in the UMI [int]
    @param sdata  Statistics of alignment between junctions and this UMI [np.array]
    """
    def __init__(self):   # this function takes no parameters as it's used by MemPool.
        self.scnt = None
        self.conf = None

        self.intv = None
        self.pe = None
        self.se = None
        self.n = 0
        self.sdata = None
        #self.sp_type = None      @param sp_type Splice type of this UMI [int]

    def destroy(self):
        pass

    def prepare(self, scnt, conf):
        self.scnt = scnt
        self.conf = conf

        self.intv = self.scnt.intv
        if conf.pair_end:
            self.pe = {}
        else:
            self.se = []

    def push_read(self, read):
        """Push one read into this count machine.
        @param read   One pysam::AlignedSegment object.
        @return       0 if success, -1 otherwise [int]
        """
        conf = self.conf
        self.n += 1
        dat = None
        if conf.pair_end:
            rname = read.query_name
            if rname in self.pe:
                dat = self.pe[rname]
            else:
                dat = self.pe[rname] = [] # TODO: fix it
        else:
            dat = self.se

        mcnt = self.scnt.mcnt
        raln = mcnt.pool_raln.get()
        if raln.align(read, mcnt) < 0:
            return(-1)
        dat.append(raln)
        return(0)

    def reset(self):
        self.intv = None
        if self.conf.pair_end:
            self.pe = {}
        else:
            self.se = []
        self.n = 0
        self.sdata = None

    # TODO: Assign different weights to UMIs with different supporting reads?
    def stat(self):
        conf = self.conf       
        n_jc = self.intv.get_n()
        rdata = np.zeros((0, n_jc, JC_R_N), dtype = ST_R_DTYPE)
        self.sdata = np.zeros((n_jc, JC_U_N), dtype = ST_U_DTYPE)

        if self.se:
            for raln in self.se:
                if raln.data:
                    r = np.zeros((1, n_jc, JC_R_N), dtype = ST_R_DTYPE)
                    for jc, aln in raln.data:
                        r[0, jc.index] = aln.sdata   # deep copy?
                    rdata = np.append(rdata, r, axis = 0)
            if rdata.size <= 0:   # this umi does not overlap with any junction.
                return(0)

            rdata = np.transpose(rdata, (1, 2, 0))  # jc x stat x read
            for i in range(rdata.shape[0]):
                sdata = rdata[i, ]
                nr_ovp = np.sum(sdata[JC_R_OVP, ])
                if nr_ovp < 1:
                    continue    # keep stat of this jc in self.sdata all 0
                self.sdata[i, JC_U_OVP] = 1

                nr_cross = np.sum(sdata[JC_R_CROSS, ])
                nr_cover = np.sum(sdata[JC_R_BOUND5] | sdata[JC_R_IN] | sdata[JC_R_BOUND3])

                if nr_cross == nr_cover:
                    self.sdata[i, JC_U_OTHER] = 1
                    continue
                nr_hit = max(nr_cross, nr_cover)
                if nr_hit < conf.min_uread or float(nr_hit) / nr_ovp < conf.min_ufrac:
                    self.sdata[i, JC_U_OTHER] = 1
                    continue                    
                if nr_cross > nr_cover:
                    self.sdata[i, JC_U_CROSS] = 1
                else:
                    self.sdata[i, JC_U_COVER] = 1
                    nr_b5 = np.sum(sdata[JC_R_BOUND5, ])  # CHECK ME! no using JC_R_BOUND5/3?
                    if nr_b5 >= conf.min_uread:
                        self.sdata[i, JC_U_BOUND5] = 1
                    nr_b3 = np.sum(sdata[JC_R_BOUND3, ])
                    if nr_b3 >= conf.min_uread:
                        self.sdata[i, JC_U_BOUND3] = 1

            return(0)

        elif self.pe:
            self.sdata = None
            return(0) #<DEV/>#pass
        else:
            self.sdata = None
            return(0)

JC_U_OVP = 0         # whether overlap or not?
JC_U_CROSS = 1
JC_U_COVER = 2
JC_U_BOUND5 = 3
JC_U_BOUND3 = 4
JC_U_OTHER = 5       # OVP and not any of [CROSS, COVER]
JC_U_N = 6

class SCount:
    """Counting for single sample
    @param mcnt     A MCount object that the SCount object belongs to.
    @param conf     A config::Config object.
    @param intv     A region::JunctionSet object.
    @param tcount   Total transcript / reads count for single sample.
    @param umi_cnt  HashMap of <str:UCount> for umi:UCount pair, mainly for 10x data.
    @param read_cnt HashMap of <str:dict> for reads_name:UCount pair, mainly for SMART-seq data.
    @param sdata    Statistics of junctions in this sample [np.array]
    """
    def __init__(self, mcnt, conf):
        self.mcnt = mcnt
        self.conf = conf

        self.intv = None   # this @p intv would be set by MCount
        self.tcount = [0] * MCNT_F_N
        self.umi_cnt = {}
        self.read_cnt = None
        if not conf.use_umi():
            self.read_cnt = mcnt.pool_ucnt.get()
        self.sdata = None
        
    def destroy(self):
        self.tcount.clear()
        self.umi_cnt.clear()  # UCount is in MCount.pool_ucnt

    def push_read(self, read):
        conf = self.conf
        ucnt = None
        if conf.use_umi():
            umi = read.get_tag(conf.umi_tag)
            if umi in self.umi_cnt:
                ucnt = self.umi_cnt[umi]
            else:
                ucnt = self.mcnt.pool_ucnt.get()
                ucnt.prepare(self, self.conf)
                self.umi_cnt[umi] = ucnt
        else:
            ucnt = self.read_cnt
        
        if self.intv.get_n() <= 0:  # when n_jc=0, just push the UMI, not the read.
            return(0)
        
        ret = ucnt.push_read(read)
        if ret < 0:
            return(-1)
        return(0)

    def reset(self):
        self.intv = None
        for i in range(len(self.tcount)):
            self.tcount[i] = 0
        self.umi_cnt.clear()
        if self.read_cnt:
            self.read_cnt.reset()

    def stat(self):
        conf = self.conf
        if not conf.use_umi():
            self.sdata = None
            return(0)  #<DEV/>#pass

        n_jc = self.intv.get_n()
        sdata = np.zeros((n_jc, JC_S_N), dtype = ST_S_DTYPE)
        for ucnt in self.umi_cnt.values():
            if ucnt.stat() < 0:
                return(-1)
            sdata[:, JC_S_OVP] += ucnt.sdata[:, JC_U_OVP]
            sdata[:, JC_S_CROSS] += ucnt.sdata[:, JC_U_CROSS]
            sdata[:, JC_S_COVER] += ucnt.sdata[:, JC_U_COVER]
            sdata[:, JC_S_OTHER] += ucnt.sdata[:, JC_U_OTHER]
        self.sdata = sdata
        return(0)

JC_S_OVP = 0
JC_S_CROSS = 1
JC_S_COVER = 2
JC_S_OTHER = 3
JC_S_N = 4

class MCount:
    """Counting for multiple samples
    @param samples    A list of sample IDs/barcodes [list of str]
    @param conf       A config::Config object.
    @param gene       A gff::Gene object.
    @param intv       A region::JunctionSet object.
    @param tcount     Array of total read / UMI counts [list of int]
    @param cell_cnt   HashMap of <str, SCount> for sample:SCount pair.
    @param jdata      Statistics of each junction [np.array]
    @param pool_ucnt  Memory Pool of UCount objects.
    @param pool_aln   Memory Pool of JCAlign objects.
    @param pool_raln  Memory Pool of JCReadAlign objects.
    @param is_reset   Has this object been reset [bool]
    """
    def __init__(self, samples, conf):
        self.samples = samples
        self.conf = conf

        self.gene = None
        self.intv = None
        self.tcount = [0] * MCNT_F_N
        self.cell_cnt = {}
        for smp in self.samples:
            if smp in self.cell_cnt:    # duplicate samples
                return(-2)
            self.cell_cnt[smp] = SCount(self, self.conf)
        self.jdata = None
        #self.sdata = None  @param sdata      Statistics of juctions in all samples [np.array]
        
        self.pool_ucnt = MemPool(UCount)
        self.pool_aln = MemPool(JCAlign)
        self.pool_raln = MemPool(JCReadAlign)

        self.is_reset = False

    def add_gene(self, gene, intv):
        conf = self.conf
        if not self.is_reset:
            self.reset()
        self.gene = gene
        self.intv = intv
        for smp, scnt in self.cell_cnt.items():
            scnt.intv = intv
        self.is_reset = False
        return(0)

    def destroy(self):
        self.intv.destroy()
        self.tcount.clear()
        for cb, scnt in self.cell_cnt.items():
            scnt.destroy()
        self.cell_cnt.clear()
        self.pool_ucnt.destroy()
        self.pool_aln.destroy()
        self.pool_raln.destroy()

    def push_read(self, read, sid = 0):
        """Push one read into this counting machine
        @param read  A pysam::AlignedSegment object.
        @param sid   Sample ID [int]
        @return      0 if success, -1 error, -2 read filtered [int]
        """
        conf = self.conf
        if conf.use_barcode():
            cb = read.get_tag(conf.cell_tag)
        elif sid < len(self.sid_list):
            cb = self.sid_list[sid]
        else:
            return(-1)

        if cb in self.cell_cnt:
            scnt = self.cell_cnt[cb]
        else:
            return(-2)

        ret = scnt.push_read(read)
        if ret < 0: 
            return(-1)
        return(0)

    def reset(self):
        if self.intv:
            self.intv.reset()
        if self.tcount:
            for i in range(len(self.tcount)):
                self.tcount[i] = 0
        if self.cell_cnt:
            for smp in self.cell_cnt:
                self.cell_cnt[smp].reset()
        if self.pool_ucnt:
            self.pool_ucnt.reset()
        if self.pool_aln:
            self.pool_aln.reset()
        if self.pool_raln:
            self.pool_raln.reset()
        self.is_reset = True

    def stat(self):
        # TODO:
        #   1. make use of info of BOUND5/BOUND3/OTHER.
        #   2. use soft cutoff instead of hard cutoff. see GATK strategy.
        #   3. include codes for denovo junctions.

        func = "MCount::stat"
        conf = self.conf
        if not conf.use_umi():
            self.jdata = None
            return(0)    #<DEV/>#pass
        
        n_cell = len(self.cell_cnt)
        n_jc = self.intv.get_n()
        if n_jc <= 0:
            self.jdata = None
            n_umi = sum([len(scnt.umi_cnt) for scnt in self.cell_cnt.values()])
            self.tcount[MCNT_F_UNSPLICED] = n_umi
            return(0)

        cdata = np.zeros((n_cell, n_jc, JC_S_N), dtype = ST_S_DTYPE)
        for i, smp in enumerate(self.samples):
            scnt = self.cell_cnt[smp]
            if scnt.stat() < 0:
                return(-1)
            cdata[i, ] = scnt.sdata
        cdata = np.transpose(cdata, (1, 0, 2))  # jc x cell x stat

        self.jdata = np.zeros((n_jc, 3), dtype = "int32")
        jstat = np.zeros((n_jc, ), dtype = ST_J_DTYPE)
        for i in range(cdata.shape[0]):
            dat = cdata[i, ]    # cell x stat
            n_jumi = np.sum(dat[:, JC_S_CROSS])
            n_jcell = np.sum(dat[:, JC_S_CROSS] > 0)
            if n_jumi >= conf.min_jumi and n_jcell >= conf.min_jcell:
                jstat[i] = 1
            self.jdata[i, 0] = n_jumi
            self.jdata[i, 1] = n_jcell
            self.jdata[i, 2] = jstat[i]

        for i, smp in enumerate(self.samples):
            scnt = self.cell_cnt[smp]
            for umi, ucnt in scnt.umi_cnt.items():
                nj_cross = np.sum(ucnt.sdata[:, JC_U_CROSS] * jstat)
                nj_cover = np.sum(ucnt.sdata[:, JC_U_COVER] * jstat)
                nj_other = np.sum(ucnt.sdata[:, JC_U_OTHER] * jstat)
                if conf.debug > 1:
                    sys.stderr.write("[D::%s] barcode=%s; umi=%s; nj_cross=%d; nj_cover=%d; nj_other=%d\n" %
                        (func, smp, umi, nj_cross, nj_cover, nj_other))
                if nj_cross > 0:
                    scnt.tcount[MCNT_F_SPLICED] += 1
                elif nj_cover > 0:
                    scnt.tcount[MCNT_F_UNSPLICED] += 1
                else:
                    scnt.tcount[MCNT_F_AMBIGUOUS] += 1
            self.tcount[MCNT_F_SPLICED] += scnt.tcount[MCNT_F_SPLICED]
            self.tcount[MCNT_F_UNSPLICED] += scnt.tcount[MCNT_F_UNSPLICED]
            self.tcount[MCNT_F_AMBIGUOUS] += scnt.tcount[MCNT_F_AMBIGUOUS]

        return(0)

ST_R_DTYPE = "int8"
ST_U_DTYPE = "int8"
ST_S_DTYPE = "int32"
ST_J_DTYPE = "int8"
ST_M_DTYPE = "int32"

# ?unannotated: reads cannot be properly aligned to any annotated transcript.
# ?ambiguous: reads can be either spliced or unspliced for one transcript, 
#     eg., pair reads are aligned to the same big exon.
MCNT_F_SPLICED = 0
MCNT_F_UNSPLICED = 1
MCNT_F_AMBIGUOUS = 2 
#MCNT_F_UNANNO = 3    
#MCNT_F_N = 4
MCNT_F_N = 3
