# count.py - feature counting utils
# Copied from cellSNP, https://github.com/single-cell-genetics/cellSNP/blob/feature-count/cellSNP/utils/count_utils.py
# Utilility functions for count reads in each feature
# Author: Yuanhua Huang
# Date: 27/04/2020

## Note, this is a very basic read counting for features

import sys
import numpy as np
import pysam

from ..config import DEBUG


# for pysam
CACHE_CHROM = None
CACHE_SAMFILE = None

def check_pysam_chrom(samFile, chrom=None):
    """Chech if samFile is a file name or pysam object, and if chrom format. 
    """
    global CACHE_CHROM
    global CACHE_SAMFILE

    if CACHE_CHROM is not None:
        if (samFile == CACHE_SAMFILE) and (chrom == CACHE_CHROM):
            return CACHE_SAMFILE, CACHE_CHROM

    if type(samFile) == str:
        ftype = samFile.split(".")[-1]
        if ftype != "bam" and ftype != "sam" and ftype != "cram" :
            print("Error: file type need suffix of bam, sam or cram.")
            sys.exit(1)
        if ftype == "cram":
            samFile = pysam.AlignmentFile(samFile, "rc")
        elif ftype == "bam":
            samFile = pysam.AlignmentFile(samFile, "rb")
        else:
            samFile = pysam.AlignmentFile(samFile, "r")

    if chrom is not None:
        if chrom not in samFile.references:
            if chrom.startswith("chr"):
                chrom = chrom.split("chr")[1]
            else:
                chrom = "chr" + chrom
        if chrom not in samFile.references:
            print("Can't find references %s in samFile" %chrom)
            return samFile, None
    
    CACHE_CHROM = chrom
    CACHE_SAMFILE = samFile
    return samFile, chrom


def fetch_reads(sam_file, region, cell_tag="CR", UMI_tag="UR", min_MAPQ=20, 
                max_FLAG=255, min_LEN=30):
    """ Fetch all reads mapped to a given region.
    Filtering is also applied, including cell and UMI tags and read mapping 
    quality.
    """
    if sam_file is None or region is None:
        if sam_file is None:
            print("Warning: samFile is None")
        if region is None:
            print("Warning: region is None")
        return np.array([]), np.array([])

    samFile, _chrom = check_pysam_chrom(sam_file, region.chrom)
    if not _chrom:
        print("Warning: chrom '%s' not in sam header" % region.chrom)
        return np.array([]), np.array([])
    
    UMIs_list, cell_list = [], []
    for _read in samFile.fetch(_chrom, region.start - 1, region.end):
        ## filtering reads
        # this might be further speed up
        overhang = sum((np.array(_read.positions) >= (region.start - 1)) *   
                       (np.array(_read.positions) <= (region.end - 1)))

        if _read.mapq < min_MAPQ or _read.flag > max_FLAG or overhang < min_LEN: 
            continue
            
        if cell_tag is not None and _read.has_tag(cell_tag) == False: 
            continue
        if UMI_tag is not None and _read.has_tag(UMI_tag) == False: 
            continue
            
        if UMI_tag is not None:
            UMIs_list.append(_read.get_tag(UMI_tag))
        if cell_tag is not None:
            cell_list.append(_read.get_tag(cell_tag))

    if len(cell_list) > 0 and len(cell_list) == len(UMIs_list):
        UMI_cell = [UMIs_list[x] + cell_list[x] for x in range(len(UMIs_list))]
        UMI_cell, idx, cnt = unique_list(UMI_cell)
        cell_list = [cell_list[x] for x in idx]
    
    cell_list_uniq, idx, read_count = unique_list(cell_list)

    return cell_list_uniq, read_count


def feature_count(sam_file, barcodes, region, reg_index, cell_tag, UMI_tag, 
    min_MAPQ, max_FLAG, min_LEN):
    """Fetch read count for a given feature.
    """
    if DEBUG:
        sys.stderr.write("[feature_count] index = %d, region = %s:%d-%d\n" % (reg_index, region.chrom, region.start, region.end))
        sys.stderr.flush()

    cell_list_uniq, read_count = fetch_reads(sam_file, region, cell_tag, UMI_tag, 
        min_MAPQ, max_FLAG, min_LEN)

    if len(cell_list_uniq) > 0:
        match_idx = id_mapping(cell_list_uniq, barcodes, uniq_ref_only=True, 
            IDs2_sorted=True)
        match_idx = np.array(match_idx, dtype = float)

        idx1 = np.where(match_idx == match_idx)[0] #remove None
        idx2 = match_idx[idx1].astype(int)
        
        out_list = []
        for j in range(len(idx2)):
            out_list.append("%d\t%d\t%d" %(reg_index, idx2[j], read_count[idx1[j]]))
        
        if DEBUG:
            sys.stderr.write("[feature_count] index = %d, len(hits) = %d\n" % (reg_index, len(idx2)))
            sys.stderr.flush()

        return "\n".join(out_list) + "\n"
    else:
        if DEBUG:
            sys.stderr.write("[feature_count] index = %d, hits = None\n" % (reg_index, ))
            sys.stderr.flush()

        return None


## The two function id_mapping() and unique_list() are copied from 
## cellSNP, https://github.com/single-cell-genetics/cellSNP/blob/purePython/cellSNP/utils/base_utils.py

def id_mapping(IDs1, IDs2, uniq_ref_only=True, IDs2_sorted=False):
    """
    Mapping IDs2 to IDs1. IDs1 (ref id) can have repeat values, but IDs2 need 
    to only contain unique ids.
    Therefore, IDs2[rv_idx] will be the same as IDs1.
    
    Parameters
    ----------
    IDs1 : array_like or list
        ids for reference.
    IDs2 : array_like or list
        ids waiting to map.
        
    Returns
    -------
    RV_idx : array_like, the same length of IDs1
        The index for IDs2 mapped to IDs1. If an id in IDs1 does not exist 
        in IDs2, then return a None for that id.
    """
    idx1 = sorted(range(len(IDs1)), key=IDs1.__getitem__)
    if IDs2_sorted:
        idx2 = range(len(IDs2))
    else:
        idx2 = sorted(range(len(IDs2)), key=IDs2.__getitem__)
    RV_idx1, RV_idx2 = [], []
    
    i, j = 0, 0
    while i < len(idx1):
        if j == len(idx2) or IDs1[idx1[i]] < IDs2[idx2[j]]:
            RV_idx1.append(idx1[i])
            RV_idx2.append(None)
            i += 1
        elif IDs1[idx1[i]] == IDs2[idx2[j]]:
            RV_idx1.append(idx1[i])
            RV_idx2.append(idx2[j])
            i += 1
            if uniq_ref_only: j += 1
        elif IDs1[idx1[i]] > IDs2[idx2[j]]:
            j += 1
            
    origin_idx = sorted(range(len(RV_idx1)), key=RV_idx1.__getitem__)
    RV_idx = [RV_idx2[i] for i in origin_idx]
    return RV_idx


def unique_list(X):
    """unique a list with index and count
    Example
    -------
    >>> unique_list([1,2,4,5,3,2,4])
    >>> ([1, 2, 3, 4, 5], [0, 1, 4, 2, 3], [1, 2, 1, 2, 1])
    """
    idx = sorted(range(len(X)), key=X.__getitem__)
    X_uniq = []
    X_count = []
    idx_uniq = []
    for i in idx:
        if len(X_uniq) == 0 or X[i] != X_uniq[-1]:
            X_uniq.append(X[i])
            X_count.append(1)
            idx_uniq.append(i)
        else:
            X_count[-1] += 1
    return X_uniq, idx_uniq, X_count
