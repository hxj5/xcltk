# Region routine
# Author: Xianjie Huang

from functools import cmp_to_key
from intervaltree import IntervalTree
import numpy as np

class Region:
    """Region class
    @param chrom   Chromosome name [str]
    @param start   1-based start pos, inclusive [int]
    @param end     1-based end pos, exclusive [int]
    @param rid     Region ID [str]
    """
    def __init__(self, chrom, start, end, rid = None):
        self.chrom = format_chrom(chrom)
        self.start = start
        self.end = end
        self._rid = rid
        self.len = self.end - self.start

    def compare(self, region):
        """Compare with another region.
        @param region   A Region object of the same class.
        @return         Comparison result represented by [int]
                          negative integer if self is smaller; 
                          0 if equal;
                          positive integer if bigger.
        """
        if self.chrom == region.chrom:
            if self.start == region.start:
                return(self.end - region.end)
            else:
                return(self.start - region.start)
        elif self.chrom < region.chrom:
            return(-1)
        else:
            return(1)

    def get_id(self):
        """Format the region id [str]
        """
        if self._rid is None:
            self._rid = "%s_%d_%d" % (self.chrom, self.start, self.end)
        return self._rid

    def get_len(self):
        return self.len

class RegionSet:
    """Region set with payload
    @param creg    Region list for each chromosome [dict]
    @param is_sorted If the regions have been sorted for each chromosome [dict]
    @param ctree   Intervaltree for each chromosome [dict]
    @param cid     Region id map [dict]
    @param n       Total number of regions [int]
    @param is_uniq Whether payloads of duplicate region ids should be discarded [bool]
    """
    def __init__(self, is_uniq = False):
        self.is_uniq = is_uniq

        self.creg = {}
        self.is_sorted = {}

        self.ctree = {}
        self.cid = {}
        self.n = 0

    def __sort_chrom_regions(self, chrom):
        if chrom in self.creg:
            if not self.is_sorted[chrom]:
                self.creg[chrom] = self.__sort_regions(self.creg[chrom])
                self.is_sorted[chrom] = True            

    def __sort_regions(self, reg_list):
        cmp_reg = lambda r1, r2: r1.compare(r2)
        return sorted(reg_list, key = cmp_to_key(cmp_reg))
    
    def add(self, region):
        """Add a new region.
        @param region  Region-like object.
        @return        0 success, 1 discarded as duplicate, -1 error [int]
        """
        chrom = format_chrom(region.chrom)
        if chrom not in self.creg:
            self.creg[chrom] = list()
            self.ctree[chrom] = IntervalTree()

        item = None
        rid = region.get_id()
        if rid in self.cid:
            if self.is_uniq:
                return(1)
            else:
                item = region
        else:
            self.cid[rid] = list()
            item = region
        if item is not None:
            self.creg[chrom].append(item)
            self.is_sorted[chrom] = False
            self.ctree[chrom][region.start:region.end] = item
            self.cid[rid].append(item)
            self.n += 1
        return(0)

    def destroy(self):
        self.reset()

    def fetch(self, chrom, start, end):
        """Fetch overlapping regions.
        @param chrom  Chromosome name [str]
        @param start  1-based start pos, inclusive [int]
        @param end    1-based end pos, exclusive [int]
        @return       All overlapping regions (not sorted) [list]
        """
        chrom = format_chrom(chrom)
        if chrom not in self.ctree:
            return([])
        tree = self.ctree[chrom]
        hits = [region for begin, end, region in tree[start:end]]
        return(hits)

    def get_n(self):
        return(self.n)

    def get_regions(self, chrom = None, sort = False):
        """Get regions for chromosome(s).
        @param chrom  Chromosome name; All chromosomes if None [str]
        @param sort   Whether to sort the regions [bool]
        @return       A list of regions [list]
        """
        ch_list = []
        if chrom is None:
            ch_list = self.creg.keys()
        else:
            chrom = format_chrom(chrom)
            if chrom in self.creg:
                ch_list = [chrom]
            else:
                return([])

        lst = []        
        for ch in ch_list:
            if sort:
                self.__sort_chrom_regions(ch)
            lst.extend(self.creg[ch])
        return(lst)

    def merge(self, rs):
        """Merge another region set
        @param rs  RegionSet object.
        @return    Num of regions merged if success, -1 if error [int]
        """
        k = 0
        reg_list = rs.get_regions()
        for region in reg_list:
            ret = self.add(region)
            if ret != 0:
                if ret < 0:
                    return(-1)
            else:
                k += 1
        return(k)

    def query(self, rid):
        """Query region(s) given its ID.
        @param rid  Region ID [str].
        @return     A list of regions; empty list if region is not in the set.
        """
        if rid in self.cid:
            return(self.cid[rid])
        return([])

    def reset(self):
        for chrom, reg_list in self.creg.items():
            reg_list.clear()
        self.creg.clear()
        self.n = 0
        self.is_sorted.clear()
        for chrom, tree in self.ctree.items():
            tree.clear()
        self.ctree.clear()
        for chrom, id_set in self.cid.items():
            id_set.clear()
        self.cid.clear()

    def sort(self):
        for chrom in self.citem:
            if not self.is_sorted[chrom]:
                self.citem[chrom] = self.__sort_items(self.citem[chrom])
                self.is_sorted[chrom] = True

class Junction(Region):
    """Splice Junction class for one gene
    @param chrom    Chromosome name [str]
    @param start    1-based start pos, inclusive [int]
    @param end      1-based end pos, exclusive [int]
    @param junc_id  ID of the junction; None to use default setting [str]
    @param tran_id  ID of the transcript the junction belongs to [str]
    @param gene_id  ID of the gene the junction belongs to [str]
    @param is_anno  Is this junction annotated [bool]
    @param exon5    A tuple of start and end pos of the 5'-adjacent exon; both pos
                      are 1-based; 5'-inclusive; 3'-exclusive [int]
    @param exon3    A tuple of start and end pos of the 3'-adjacent exon; both pos
                      are 1-based; 5'-inclusive; 3'-exclusive [int]
    #"""
    def __init__(self, chrom, start, end, junc_id,
                 tran_id, gene_id,
                 is_anno,
                 exon5 = None, exon3 = None):
        super().__init__(chrom, start, end, junc_id)
        self.tran_id = tran_id
        self.gene_id = gene_id
        self.is_anno = is_anno
        self.exon5 = exon5      # exon5 and exon3 for #DEV# version.
        self.exon3 = exon3

        self.index = None     # NOTE that the index must start from 0. see mcount::MCount.

    #TODO: could be extended to include transcript (id_map) info.
    #def get_id(self):
    #    pass
    
class JunctionSet(RegionSet):
    """Splice Junction set for one gene"""
    def __init__(self, is_uniq = False):
        super().__init__(is_uniq)

    def add(self, junction):
        ret = super().add(junction)
        if ret != 0:
            return(ret)
        else:
            junction.index = self.n - 1
            return(0)

    def get_junctions(self, sort = True):
        return self.get_regions(chrom = None, sort = sort)

    def get_len_quantile(self, quantile, anno_only = True):
        """Get length quantile of the regions
        @param quantile   [float]
        @param anno_only  Whether to use annotated regions only [bool]
        @return The length quantile or None [float]
        """
        regions = self.get_junctions()
        if anno_only:
            reg_len = [reg.get_len() for reg in regions if reg.is_anno]
        else:
            reg_len = [reg.get_len() for reg in regions]
        if not reg_len:
            return None
        return np.quantile(reg_len, quantile)

class XonSet(RegionSet):
    """The exon/intron level region set"""
    def __init__(self, is_uniq = False):
        super().__init__(is_uniq)

def format_chrom(chrom):
    return chrom[3:] if chrom.lower().startswith("chr") else chrom

REG_EXON = 1
REG_INTRON = 2

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 3:
        sys.stdout.write("Usage: %s <gtf|gff> <chrom:start-end>\n" % sys.argv[0])
        sys.exit(1)
    in_fn = sys.argv[1]
    region = sys.argv[2]
    chrom = region.split(":")[0]
    start = int(region.split(":")[1].split("-")[0])
    end = int(region.split("-")[1]) + 1

    from gff import gff_load
    ret, gene_set = gff_load(in_fn)
    if ret < 0:
        sys.stderr.write("Error: gff_load failed.\n")
        sys.exit(3)
    genes = gene_set.get_genes()

    rs = RegionSet()
    for g in genes:
        reg = Region(g.chrom, g.start, g.end + 1, g.gene_id)
        rs.add(reg)
    hits = rs.fetch(chrom, start, end + 1)
    for h in hits:
        sys.stdout.write("%s,%d,%d,%s\n" % (h.chrom, h.start, h.end, h.rid))
