# gfeature.py - genomic features, supporting interval query.


from ...utils.grange import Region, RegionSet



class SNP(Region):
    """Phased SNP

    Attributes
    ----------
    chrom : str
        Chromosome name.
    pos : int
        1-based position.
    ref : str
        The ref base.
    alt : str
        The alt base.
    ref_idx : int
        The GT index for ref base, 0 or 1.
    alt_idx : int
        The GT index for alt base, 1 or 0.   
    """
    def __init__(self, chrom, pos, ref, alt, ref_idx, alt_idx):
        super().__init__(chrom, pos, pos + 1)
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.ref_idx = ref_idx
        self.alt_idx = alt_idx
        self.gt = {ref:ref_idx, alt:alt_idx}

    def get_id(self):
        return "%s_%d" % (self.chrom, self.pos)

    def get_region_allele_index(self, base):
        return self.gt[base] if base in self.gt else -1


    
class SNPSet(RegionSet):
    """A set of phased SNPs"""
    def __init__(self, is_uniq = False):
        super().__init__(is_uniq)

    def add(self, snp):
        return super().add(snp)


    
class BlockRegion(Region):
    """Block Region

    Attributes
    ----------
    chrom : str
        Chromosome name.
    start : int
        1-based start pos, inclusive.
    end : int
        1-based end pos, exclusive.
    name : str
        Name of the block.
    snp_list : list
        A list of SNPs (`SNP` objects) located within the block.
    """
    def __init__(self, chrom, start, end, name = None, snp_list = None):
        super().__init__(chrom, start, end)
        self.name = name
        self.snp_list = snp_list
