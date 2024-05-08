# mcount.py - Reads Counting
# Author: Xianjie Huang

from ...utils.sam import get_query_bases


class UCount:
    """Counting Unit

    Parameters
    ----------
    scnt : SCount object
        A SCount object that the UCount object belongs to.
    conf : config::Config object
        Configuration.
    allele : str
        The allele for the query SNP in this UMI.
    """
    def __init__(self):
        self.scnt = None
        self.conf = None
        self.allele = None

    # here we use `prepare` function to set `UCount`, instead of initialize
    # directly in `__init__` function, in case a memory pool of `UCount`
    # is used.
    def prepare(self, scnt, conf):
        self.scnt = scnt
        self.conf = conf

    def push_read(self, read):
        """Push one read into this count machine.
        
        Parameters
        ----------
        read : pysam::AlignedSegment object
            A BAM read to be counted.

        Returns
        -------
        int
            0 if success, -1 otherwise.
        """
        snp = self.scnt.mcnt.snp
        try:
            idx = read.positions.index(snp.pos - 1)
        except:
            self.allele = None
        else:
            bases = get_query_bases(read, full_length = False)
            self.allele = bases[idx].upper()
        return(0)

    def stat(self):
        return(0)


class SCount:
    """Counting for single sample

    Parameters
    ----------
    mcnt : MCount object
        A MCount object that the SCount object belongs to.
    conf : config::Config object
        Configuration.
    tcount : list
        Total read / UMI counts for A/C/G/T/N bases, only for this 
        sample [list of int; 5 elements].
    umi_cnt : dict
        HashMap of <str:UCount> for umi:UCount pair, mainly for 10x data.
    """
    def __init__(self, mcnt, conf):
        self.mcnt = mcnt
        self.conf = conf

        self.tcount = [0] * 5
        self.umi_cnt = {}

    def push_read(self, read):
        conf = self.conf
        umi = None
        if conf.use_umi():
            umi = read.get_tag(conf.umi_tag)
        else:
            umi = read.query_name
        if not umi:
            return(0)
        if umi in self.umi_cnt:
            return(0)
        else:
            ucnt = UCount()
            ucnt.prepare(self, conf)
            self.umi_cnt[umi] = ucnt
            ret = ucnt.push_read(read)
            if ret < 0:
                return(-1)
            return(0)

    def reset(self):
        for i in range(len(self.tcount)):
            self.tcount[i] = 0
        self.umi_cnt.clear()

    def stat(self):
        for ucnt in self.umi_cnt.values():
            if ucnt.stat() < 0:
                return(-1)
            allele = ucnt.allele
            if allele:
                idx = self.mcnt.base_idx["N"] 
                if allele in self.mcnt.base_idx:
                    idx = self.mcnt.base_idx[allele]
                self.tcount[idx] += 1
        return(0)


class MCount:
    """Counting for multiple samples

    Parameters
    ----------
    samples : list
        A list of cell barcodes or sample IDs [list of str].
    conf : config::Config object
        Configuration
    snp : region::SNP object
        The SNP beging pileuped.
    tcount : list
        Total read / UMI counts for A/C/G/T/N bases, aggregated for all 
        samples [list of int; 5 elements].
    base_idx : dict
        The mapping from base (str) to index (int) for `tcount`.
    cell_cnt : dict
        HashMap of <str, SCount> for sample:SCount pair.
    is_reset : boot
        Has this object been reset.
    """
    def __init__(self, samples, conf):
        self.samples = samples
        self.conf = conf

        self.snp = None
        self.tcount = [0] * 5
        self.base_idx = {"A":0, "C":1, "G":2, "T":3, "N":4}
        self.cell_cnt = {}
        for smp in self.samples:
            if smp in self.cell_cnt:    # duplicate samples
                return(-2)
            self.cell_cnt[smp] = SCount(self, self.conf)
        self.is_reset = False

    def add_snp(self, snp):
        if not self.is_reset:
            self.reset()
        self.snp = snp
        self.is_reset = False
        return(0)

    def push_read(self, read, sid = None):
        """Push one read into this counting machine.

        Parameters
        ----------
        read : pysam::AlignedSegment object
            A BAM read to be counted.
        sid : str
            The ID of the sample that the read belongs to. 
            Set to `None` if cell barcodes are used.

        Returns
        -------
        int
            0 if success, -1 error, -2 read filtered.
        """
        conf = self.conf
        if conf.use_barcodes():
            smp = read.get_tag(conf.cell_tag)
        else:
            smp = sid
        scnt = None
        if smp in self.cell_cnt:
            scnt = self.cell_cnt[smp]
        else:
            return(-2)

        ret = scnt.push_read(read)
        if ret < 0: 
            return(-1)
        return(0)

    def reset(self):
        if self.tcount:
            for i in range(len(self.tcount)):
                self.tcount[i] = 0
        if self.cell_cnt:
            for smp in self.cell_cnt:
                self.cell_cnt[smp].reset()
        self.is_reset = True

    def stat(self):
        for smp in self.samples:
            scnt = self.cell_cnt[smp]
            if scnt.stat() < 0:
                return(-1)
            for j in range(len(self.tcount)):
                self.tcount[j] += scnt.tcount[j]
        return(0)