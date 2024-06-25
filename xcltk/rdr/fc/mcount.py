# mcount.py - counting machine for features.


class SCount:
    """Counting for single sample

    Attributes
    ----------
    mcnt : MCount object
        A MCount object that the SCount object belongs to.
    conf : config::Config object
        Configuration.
    tcount : int
        Total read / UMI counts, only for this sample.
    umi_set : set
        The set of UMIs belonging to the feature and the sample.
    is_rest : bool
        Has this object been reset?
    """
    def __init__(self, mcnt, conf):
        self.mcnt = mcnt
        self.conf = conf

        self.tcount = 0
        self.umi_set = set()
        self.is_reset = False

    def mark_reset_false(self):
        self.is_reset = False

    def mark_reset_true(self):
        self.is_reset = True

    def push_read(self, read):
        conf = self.conf
        umi = None
        if conf.use_umi():
            umi = read.get_tag(conf.umi_tag)
        else:
            umi = read.query_name
        if umi:
            self.umi_set.add(umi)
        return(0)

    def reset(self):
        if self.is_reset:
            return
        self.tcount = 0
        self.umi_set.clear()
        self.mark_reset_true()

    def stat(self):
        self.tcount = len(self.umi_set)
        return(0)


class MCount:
    """Counting for multiple samples

    Attributes
    ----------
    samples : list
        A list of cell barcodes or sample IDs [list of str].
    conf : config::Config object
        Configuration
    tcount : list
        Total read / UMI counts, aggregated for all samples.
    cell_cnt : dict
        HashMap of <str, SCount> for sample:SCount pair.
    is_reset : boot
        Has this object been reset.
    """
    def __init__(self, samples, conf):
        self.samples = samples
        self.conf = conf

        self.tcount = 0
        self.cell_cnt = {}
        for smp in self.samples:
            if smp in self.cell_cnt:    # duplicate samples
                return(-2)
            self.cell_cnt[smp] = SCount(self, self.conf)
        self.is_reset = False

    def add_region(self, reg):
        self.reset()
        self.mark_reset_false()
        return(0)
    
    def mark_reset_false(self, recursive = True):
        self.is_reset = False
        if recursive:
            for scnt in self.cell_cnt.values():
                scnt.mark_reset_false()

    def mark_reset_true(self, recursive = False):
        self.is_reset = True
        if recursive:
            for scnt in self.cell_cnt.values():
                scnt.mark_reset_true()

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
        if self.is_reset:
            return
        self.tcount = 0
        if self.cell_cnt:
            for scnt in self.cell_cnt.values():
                scnt.reset()
        self.mark_reset_true()

    def stat(self):
        for smp in self.samples:
            scnt = self.cell_cnt[smp]
            if scnt.stat() < 0:
                return(-1)
            self.tcount += scnt.tcount
        return(0)