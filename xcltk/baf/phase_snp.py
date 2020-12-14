#!/usr/bin/env python
#-*-coding:utf8-*-
#this script is aimed to aggregate SNPs into haplotype blocks to get BAF matrices within each cell.
#Author: Xianjie Huang 

import os
import sys
import getopt
import gzip
from ..utils.base import assert_path_exists, log
from .config import APP

def __format_chrom(chrom):
    """
    @abstract    Format chrom name to keep the chroms used in this script in the same style
    @param chrom Chrom name [str]
    @return      Formatted chrom name [str]
    """
    return chrom[3:] if chrom.startswith("chr") else chrom

def __load_snp_mtx(fn):
    """
    @abstract  Load data from SNP AD/DP mtx
    @param fn  Path to mtx file [str]
    @return    A tuple of four elements if success, None otherwise [tuple]:
                 - number of snps [int]
                 - number of cells [int]
                 - number of records [int]
                 - a dict of {snp_idx:{cell_idx:depth, }} pairs [dict]
    """
    if not fn or not os.path.isfile(fn):
        return None
    cols = []
    if fn.endswith(".gz"):
        cols = [line[:-1].split("\t")[:3] for line in gzip.open(fn, "rt")]
    else:
        cols = [line[:-1].split("\t")[:3] for line in open(fn, "r")]
    if len(cols) < 3 or len(cols[2]) < 3:
        return None
    nsnp, ncell, nrecord = [int(i) for i in cols[2][:3]]
    snp_cell = {}
    for c in cols[3:]:
        assert len(c) >= 3, "too few columns in snp mtx"
        depth = int(c[2])
        if depth > 0:
            snp_cell.setdefault(c[0], {})[c[1]] = depth
    return (nsnp, ncell, nrecord, snp_cell)

def __load_phase(fn):
    """
    @abstract  Load data from phase file
    @param fn  Path to phase file [str]
    @return    A tuple of three elements if success, None otherwise [tuple]:
                 - number of total snps [int]
                 - number of valid snps whose one allele is 0 and the other is 1 [int]
                 - a dict of {chrom:[(pos, allele1, allele2, snp_idx),]} pairs [dict]
    """
    if not fn or not os.path.isfile(fn):
        return None
    cols = []
    if fn.endswith(".gz"):
        cols = [line[:-1].split("\t")[:4] for line in gzip.open(fn, "rt")]
    else:
        cols = [line[:-1].split("\t")[:4] for line in open(fn, "rt")]
    phases = {}
    i, j = 0, 0
    for c in cols:
        assert len(c) >= 4, "too few columns in phase file"
        i += 1
        if (c[2] == "0" and c[3] == "1") or (c[2] == "1" and c[3] == "0"):
            chrom = __format_chrom(c[0])
            phases.setdefault(chrom, []).append((int(c[1]), c[2], c[3], str(i)))
            j += 1
    return (i, j, phases)

def __load_region(fn):
    """
    @abstract  Load data from region file.
    @param fn  Path to region file [str]
    @return    A tuple of two elements if success, None otherwise [tuple]:
                 - number of blocks [int]
                 - a dict of {chrom:[(start, end, reg_idx),]} pairs [dict]
    """
    if not fn or not os.path.isfile(fn):
        return None
    cols = []
    if fn.endswith(".gz"):
        cols = [line[:-1].split("\t")[:3] for line in gzip.open(fn, "rt")]
    else:
        cols = [line[:-1].split("\t")[:3] for line in open(fn, "rt")]
    regions = {}
    for i, c in enumerate(cols):    # here enumerate is efficient as region file is usually small.
        assert len(c) >= 3, "too few columns in region file"
        chrom = __format_chrom(c[0])
        regions.setdefault(chrom, []).append((int(c[1]), int(c[2]), i + 1))
    return (i + 1, regions)

def __get_block_cell(snp_ad, snp_dp, phase, blocks):
    """
    @abstract        Get block-cell AD & DP matrices.
    @param snp_ad    SNP AD mtx, A dict of {snp_idx:{cell_idx:depth, }} pairs [dict]
    @param snp_dp    SNP DP mtx, A dict of {snp_idx:{cell_idx:depth, }} pairs [dict]
    @param phase     A dict of {chrom:[(pos, allele1, allele2, snp_idx),]} pairs,
                     every array of chr should be sorted by pos already [dict]
    @param blocks    A dict of {chrom:[(start, end, reg_idx),]} pairs,
                     every array of chr should be sorted by start pos already [dict]
    @return          A dict of {reg_idx:{cell_idx:{ad:ad_depth, dp:dp_depth}, }} pairs if success, None otherwise [dict]
    @note            1. SNP AD & DP mtx should have been checked to make sure each snp-cell
                        record whose depth > 0 in AD should also exist in DP.
                     2. The two alleles of each snp should be one is 0 and the other is 1.
    """
    if not (snp_ad and snp_dp and phase and blocks):
        return None
    block_cell = {}
    for chrom, reg_dat in blocks.items():
        ph_dat = phase.get(chrom, [])
        if not ph_dat:
            continue
        ph_idx = 0
        nph = len(ph_dat)
        for r in reg_dat:
            start, end, reg_idx = r[:3]
            _idx = ph_idx - 1
            while _idx >= 0:
                pos, allele1, allele2, snp_idx = ph_dat[_idx][:4]
                if pos < start:
                    break
                elif pos <= end:
                    _idx -= 1
                    dp_dat = snp_dp.get(snp_idx, {})
                    if not dp_dat:
                        continue
                    ad_dat = snp_ad.get(snp_idx, {})
                    for cell_idx, dp_depth in dp_dat.items():
                        ad_depth = ad_dat.get(cell_idx, 0)
                        block_cell.setdefault(reg_idx, {}).setdefault(cell_idx, {"ad":0, "dp":0})
                        block_cell[reg_idx][cell_idx]["ad"] += ad_depth if allele1 == "1" else dp_depth - ad_depth
                        block_cell[reg_idx][cell_idx]["dp"] += dp_depth
                else:
                    _idx -= 1
            while ph_idx < nph:    # donot use range() here to save memory
                pos, allele1, allele2, snp_idx = ph_dat[ph_idx][:4]
                if pos > end:
                    break
                elif pos < start: 
                    ph_idx += 1
                    continue
                else:          # this snp belongs to the block
                    ph_idx += 1            
                    dp_dat = snp_dp.get(snp_idx, {})
                    if not dp_dat:
                        continue
                    ad_dat = snp_ad.get(snp_idx, {})
                    for cell_idx, dp_depth in dp_dat.items():
                        ad_depth = ad_dat.get(cell_idx, 0)
                        block_cell.setdefault(reg_idx, {}).setdefault(cell_idx, {"ad":0, "dp":0})
                        block_cell[reg_idx][cell_idx]["ad"] += ad_depth if allele1 == "1" else dp_depth - ad_depth
                        block_cell[reg_idx][cell_idx]["dp"] += dp_depth
    return block_cell
                    
def __output_block_mtx(block_cell, nblock, ncell, block_ad_file, block_dp_file, _gzip = 0):
    """
    @abstract            Output block-cell AD & DP matrices to mtx file
    @param block_cell    A dict of {reg_idx:{cell_idx:{ad:ad_depth, dp:dp_depth}, }} pairs [dict]
    @param nblock        Number of total blocks [int]
    @param ncell         Number of total cells [int]
    @param block_ad_file Path to block AD mtx file [str]
    @param block_dp_file Path to block DP mtx file [str]
    @param _gzip         If the output files need to be gziped: 0, no; 1, yes [int]
    @return              0 if success, -1 otherwise [int]
    """
    if not (block_cell and block_ad_file and block_dp_file):
        return -1

    # count number of total records in block AD & DP matrices
    nrec_ad, nrec_dp = 0, 0
    for reg_idx, reg_dat in block_cell.items():
        for cell_idx, cell_dat in reg_dat.items():
            if cell_dat["ad"] > 0: nrec_ad += 1
            if cell_dat["dp"] > 0: nrec_dp += 1

    # output mtx header
    ad_fp = gzip.open(block_ad_file, "wb") if _gzip else open(block_ad_file, "w")
    dp_fp = gzip.open(block_dp_file, "wb") if _gzip else open(block_dp_file, "w")
    header = "%%MatrixMarket matrix coordinate integer general\n%"
    ad_fp.write("%s\n%d\t%d\t%d\n" % (header, nblock, ncell, nrec_ad))
    dp_fp.write("%s\n%d\t%d\t%d\n" % (header, nblock, ncell, nrec_dp))

    # output mtx records
    sorted_reg_idx = sorted(block_cell.keys(), key = lambda x: int(x))
    for reg_idx in sorted_reg_idx:
        reg_dat = block_cell[reg_idx]
        sorted_cell_idx = sorted(reg_dat.keys(), key = lambda x: int(x))
        for cell_idx in sorted_cell_idx:
            cell_dat = reg_dat[cell_idx]
            if cell_dat["ad"] > 0:
                ad_fp.write("%s\t%s\t%d\n" % (reg_idx, cell_idx, cell_dat["ad"]))
            if cell_dat["dp"] > 0:
                dp_fp.write("%s\t%s\t%d\n" % (reg_idx, cell_idx, cell_dat["dp"]))
    ad_fp.close()
    dp_fp.close()
    return 0

def __phase_snp2block(sid, snp_ad_file, snp_dp_file, phase_file, region_file, out_dir):
    """
    @abstract           Phase (aggregate) SNP into haplotype blocks, output AD & DP mtx of each block.
    @param sid          Sample ID [str]
    @param snp_ad_file  Path to SNP AD file [str]
    @param snp_dp_file  Path to SNP DP file [str]
    @param phase_file   Path to the phase file, 4 columns: <chr> <pos> <allele1:0|1> <allele2:0|1>;
                        pos is 1-based and the <chr> column is of the same format with region file [str]
    @param region_file  Path to region file, 3 columns: <chr> <start> <end>; 
                        Both start and stop are 1-based and included [str]
    @param out_dir      Path to output dir [str]
    @return             Void
    """
    # check args
    if not sid:
        raise ValueError("Sample ID needed.")
    assert_path_exists(snp_ad_file, "file", "SNP AD mtx") 
    assert_path_exists(snp_ad_file, "file", "SNP DP mtx")
    assert_path_exists(region_file, "file", "Region file")
    if not out_dir:
        raise ValueError("Out dir needed.")
    elif not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    # load SNP AD & DP mtx
    log("loading SNP AD mtx ...")
    res_ad = __load_snp_mtx(snp_ad_file)
    assert res_ad, "failed to load SNP AD mtx."
    nsnp_ad, ncell_ad, nrec_ad, snp_ad = res_ad[:4]
    log("AD mtx header: nsnp = %d, ncell = %d, nrec = %d" % (nsnp_ad, ncell_ad, nrec_ad))
    log("AD mtx: #uniq SNPs = %d" % len(snp_ad))

    log("loading SNP DP mtx ...")
    res_dp = __load_snp_mtx(snp_dp_file)
    assert res_dp, "failed to load SNP DP mtx."
    nsnp_dp, ncell_dp, nrec_dp, snp_dp = res_dp[:4]
    log("DP mtx header: nsnp = %d, ncell = %d, nrec = %d" % (nsnp_dp, ncell_dp, nrec_dp))
    log("DP mtx: #uniq SNPs = %d" % len(snp_dp))

    # check AD & DP mtx
    # all snp-cell combinations of AD should also exist in DP
    log("check AD & DP mtx ...")
    assert (nsnp_ad == nsnp_dp and ncell_ad == ncell_dp and nrec_ad <= nrec_dp), \
           "invalid AD or DP mtx header"
    valid = True
    for snp_idx, ad_dat in snp_ad.items():
        if snp_idx not in snp_dp:
            valid = False
            break
        for cell_idx, depth in ad_dat.items():
            if cell_idx not in snp_dp[snp_idx]:
                valid = False
                break
    assert valid, "AD mtx has snp-cell record not exist in DP mtx"
    ncell = ncell_dp

    # load phase file
    log("loading phase file ...")
    res_phase = __load_phase(phase_file)
    assert res_phase, "failed to load phase file."
    nsnp_phase, nsnp_phase_valid, phase = res_phase[:3]
    log("Phase file: total SNPs = %d, valid SNPs = %d" % (nsnp_phase, nsnp_phase_valid))
    assert nsnp_phase == nsnp_ad, "#snp of phase file not equal to the value in AD mtx header"
    assert nsnp_phase_valid >= len(snp_dp), "some snp_idx of DP not exist in phase file"

    # remove SNPs, of phase file, whose snp_idx does not exist in DP mtx.
    if nsnp_phase_valid > len(snp_dp):
        log("remove SNPs, of phase file, whose snp_idx does not exist in DP mtx ...")
        n, _ph = 0, {}
        for chrom, data in phase.items():
            _ph[chrom] = [d for d in data if d[3] in snp_dp]
            n += len(_ph[chrom])
        phase = _ph
        log("Phase file: total SNPs after removing those whose snp_idx not in DP = %d" % n)

    # sort snps by pos for phase file
    log("sort snps by pos for phase file ...")
    for chrom, data in phase.items():
        data.sort(key = lambda x: x[0])

    # load region file
    log("loading region file ...")
    res_block = __load_region(region_file)
    assert res_block, "failed to load region file."
    nblock, blocks = res_block[:2]

    # sort snps and blocks by pos first for region file
    log("sort blocks by start pos for region file ...")
    for chrom, data in blocks.items():
        data.sort(key = lambda x: x[0])
        # check if regions are non-overlapping
        #i, n = 0, len(data)
        #while i < n - 1:
        #    assert data[i][1] < data[i + 1][0], "some regions are overlapping"
        #    i += 1

    # get block-cell AD & DP matrices
    log("get block-cell AD & DP matrices ...")
    block_cell = __get_block_cell(snp_ad, snp_dp, phase, blocks)
    assert block_cell, "failed to get block-cell AD & DP matrices."

    # output block AD & DP matrices
    # you can use gzip 
    log("output block AD & DP matrices ...")
    block_ad_file = os.path.join(out_dir, sid + ".block.AD.mtx")
    block_dp_file = os.path.join(out_dir, sid + ".block.DP.mtx")
    ret = __output_block_mtx(block_cell, nblock, ncell, block_ad_file, block_dp_file, 0)
    assert ret == 0, "failed to output block AD & DP matrices."

    log("All Done")

def __usage(fp = sys.stderr):
    msg =  "\n"
    msg += "Usage: %s %s [options]\n" % (APP, COMMAND)
    msg += "\n"                                                  \
           "Options:\n"                                           \
           "  --sid STR       Sample ID.\n"                        \
           "  --snpAD FILE    Path to the SNP AD mtx, snp_idx and cell_idx are both 1-based.\n"   \
           "  --snpDP FILE    Path to the SNP DP mtx, snp_idx and cell_idx are both 1-based.\n"   \
           "  --phase FILE    Path to the SNP phase file, 4 columns:\n"                           \
           "                  <chr> <pos> <allele1:0|1> <allele2:0|1>; pos is 1-based.\n"         \
           "  --region FILE   Path to region file, 3 columns: <chr> <start> <end>;\n"             \
           "                  Both start and stop are 1-based and included.\n"                    \
           "  --outdir DIR    Path to output dir.\n"                                             \
           "  -h, --help      Print this message.\n"                                             \
           "\n"
    fp.write(msg)

def phase_snp(argv):
    # parse and check command line
    if len(argv) < 3:
        __usage(sys.stderr)
        sys.exit(1)
        
    opts, args = getopt.getopt(argv[2:], "-h", ["help", "sid=", "snpAD=", "snpDP=", "phase=", "region=", "outdir="])
    sid = snp_ad_file = snp_dp_file = phase_file = region_file = out_dir = None
    for op, val in opts:
        if op in ("--sid", ): sid = val
        elif op in ("--snpAD", ): snp_ad_file = val
        elif op in ("--snpDP", ): snp_dp_file = val
        elif op in ("--phase", ): phase_file = val
        elif op in ("--region", ): region_file = val
        elif op in ("--outdir", ): out_dir = val
        elif op in ("-h", "--help"): __usage(sys.stderr); sys.exit(1)
        else: sys.stderr.write("invalid option: %s\n" % op); sys.exit(1)

    __phase_snp2block(sid, snp_ad_file, snp_dp_file, phase_file, region_file, out_dir)

COMMAND = "phase_snp"

if __name__ == "__main__":
    phase_snp(sys.argv)

