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

class SNP(Region):
    """Phased SNP
    @param chrom    Chromosome name [str]
    @param pos      1-based position [int]
    @param ref      The ref base [str]
    @param alt      The alt base [str]
    @param ref_idx  The GT index for ref base, 0 or 1 [int]
    @param alt_idx  The GT index for alt base, 1 or 0 [int]   
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
    @param chrom    Chromosome name [str]
    @param start    1-based start pos, inclusive [int]
    @param end      1-based end pos, exclusive [int]
    @param name     Name of the block [str]
    @param snp_list A list of SNPs located within the block [list of SNP objects]
    """
    def __init__(self, chrom, start, end, name = None, snp_list = None):
        super().__init__(chrom, start, end)
        self.name = name
        self.snp_list = snp_list

def format_chrom(chrom):
    return chrom[3:] if chrom.lower().startswith("chr") else chrom

REG_EXON = 1
REG_INTRON = 2

if __name__ == "__main__":
    pass
