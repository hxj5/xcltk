# grange.py - genomic range/interval routine

from functools import cmp_to_key
from intervaltree import IntervalTree
import numpy as np


class Region:
    """Region class

    Attrbutes
    ---------
    chrom : str
        Chromosome name.
    start : int
        1-based start pos, inclusive.
    end : int
        1-based end pos, exclusive.
    rid : str
        Region ID.
    """
    def __init__(self, chrom, start, end, rid = None):
        self.chrom = format_chrom(chrom)
        self.start = start
        self.end = end
        self._rid = rid
        self.len = self.end - self.start

    def compare(self, region):
        """Compare with another region.

        Parameters
        ----------
        region : `Region`
            The region to be compared with `self`.
        
        Returns
        -------
        int
            Comparison result
            - negative integer if `self` is smaller; 
            - 0 if equal;
            - positive integer if bigger.
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

    Attributes
    ----------
    creg : dict
        Region list for each chromosome.
    is_sorted : bool
        Whether the regions have been sorted for each chromosome.
    ctree : dict
        Intervaltree for each chromosome.
    cid : dict
        Region ID map.
    n : int
        Total number of regions.
    is_uniq : bool
        Whether payloads of duplicate region ids should be discarded.
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

        Parameters
        ----------
        region : `Region` object.
            The region to be added.
        
        Returns
        -------
        int
            return code. 0 success, 1 discarded as duplicate, -1 error.
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

        Parameters
        ----------
        chrom : str
            Chromosome name.
        start : int
            1-based start pos, inclusive.
        end : int
            1-based end pos, exclusive.

        Returns
        -------
        list
            All overlapping regions (not sorted).
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

        Parameters
        ----------
        chrom : str
            Chromosome name; set to `None` to use all chromosomes.
        sort : bool
            Whether to sort the regions.

        Returns
        -------
        list
            A list of regions.
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

        Parameters
        ----------
        rs : RegionSet object
            The set of regions to be merged.
        
        Returns
        -------
        int
            Number of regions merged if success, -1 if error.
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

        Parameters
        ----------        
        rid : str
            Region ID.

        Returns
        -------
        list            
            A list of regions; empty list if region is not in the set.
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


def format_chrom(chrom):
    return chrom[3:] if chrom.lower().startswith("chr") else chrom


if __name__ == "__main__":
    pass
