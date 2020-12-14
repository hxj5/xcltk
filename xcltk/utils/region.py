# Utils for Genome Region Processing
# Author: Xianjie Huang 

import os
import sys
import gzip 
from .gtf import load_genes

CHROM_LEN_HG19 = (249250621, 243199373, 198022430, 191154276, 180915260, 171115067,  # chr1 - chr6
                  159138663, 146364022, 141213431, 135534747, 135006516, 133851895,  # chr7 - chr12
                  115169878, 107349540, 102531392, 90354753, 81195210, 78077248,     # chr13 - chr 18
                  59128983, 63025520, 48129895, 51304566, 155270560, 59373566)       # chr19 - chrY

CHROM_LEN_HG38 = (248956422, 242193529, 198295559, 190214555, 181538259, 170805979,  # chr1 - chr6
                  159345973, 145138636, 138394717, 133797422, 135086622, 133275309,  # chr7 - chr12
                  114364328, 107043718, 101991189, 90338345, 83257441, 80373285,     # chr13 - chr18
                  58617616, 64444167, 46709983, 50818468, 156040895, 57227415)       # chr19 - chrY

class Region:
    '''
    @abstract      A wrapper for different region formats, e.g. bed/gff
    @param chrom   Chromosome name [str]
    @param start   Start pos of this region: 1-based, included [int]
    @param stop    End pos of this region: 1-based, included [int]
    @param _id     Region ID [str]
    '''
    def __init__(self, chrom=None, start=None, stop=None, _id=None):
        self.chrom = chrom
        self.start = start
        self.stop = stop
        self.id = _id
    
def load_regions(reg_file, reg_type):
    """
    @abstract        Get list of regions (1-based) from input file.
    @param reg_file  Path to input region file [str]
    @param reg_type  Region type, one of bed|gff|tsv [str]
    @return          A list of Region objects if success, None otherwise [list]
    """
    if not reg_file or not os.path.isfile(reg_file) or not reg_type:
        return None
    reg_type = reg_type.lower()
    if reg_type == "gff":
        genes = load_genes(reg_file, tranTag="", exonTag="")
        return [Region(g.chrom, g.start, g.stop, g.geneID) for g in genes]
    elif reg_type not in ("bed", "tsv"):
        return None

    lines = None
    if os.path.splitext(reg_file)[1] in (".gz", ".gzip"):
        with gzip.open(reg_file, "rt") as zfp:
            lines = zfp.readlines()
    else:
        with open(reg_file, "r") as fp:
            lines = fp.readlines()
    
    reg_list = []
    for idx, ln in enumerate(lines):   # it's probably fine to use enumerate here for bed files are usually small.
        parts = ln[:-1].split("\t")
        try:
            start, stop = int(parts[1]), int(parts[2])
        except (IndexError, ValueError) as e:
            print("Error: invalid bed record in No.%d line: %s" % (idx + 1, str(e)))
            return None
        if reg_type == "bed":
            start += 1
        _id = "%s:%d-%d" % (parts[0], start, stop)
        reg_list.append(Region(parts[0], start, stop, _id))

    return reg_list

def bed2reg(bed_file):
    """
    @abstract        Get list of regions (1-based) from bed file
    @param bed_file  Path to bed file [str]
    @return          A list of Region objects if success, None otherwise [list]
    """
    return load_regions(bed_file, "bed")

def gff2reg(gff_file):
    """
    @abstract        Load regions from gff file.
    @param gff_file  Path to gff file [str]
    @return          A list of Region objects if success, None otherwise [list]
    """
    return load_regions(bed_file, "gff")

def tsv2reg(tsv_file):
    """
    @abstract        Get list of regions (1-based) from tsv file
    @param tsv_file  Path to tsv file [str]
    @return          A list of Region objects if success, None otherwise [list]
    """
    return load_regions(tsv_file, "tsv")

def chr2reg(chrom_name, chrom_len, bin_size):
    """
    @abstract         Split the chromosome to several fixed-size bins.
    @param chrom_name Chromosome name [str]
    @param chrom_len  Length of the chrom [int]
    @param bin_size   Size of the bin in bp [int]
    @return           A list of Region objects if success, None otherwise [list]
    """
    try:
        chrom_len, bin_size = int(chrom_len), int(bin_size)
    except:
        return None
    if chrom_len <= 0 or bin_size <= 0: 
        return None
    nbins = chrom_len // bin_size
    if chrom_len % bin_size != 0: nbins += 1

    reg_list = []
    for i in range(nbins):
        start = bin_size * i + 1
        stop = bin_size * (i + 1)
        _id = "%s:%d-%d" % (chrom_name, start, stop)
        reg_list.append(Region(chrom_name, start, stop, _id))
    
    return reg_list

def get_fixsize_reg_from_input_len(chroms, bin_size):
    """
    @abstract        Create fixed-size bins for certain chroms acoording to pre-defined chrom length. 
    @param chroms    Dict with <chrom_name:chrom_len> pairs [dict]
    @param bin_size  Size of bin in kb [int]
    @return          A list of Region objects if success, None otherwise [list]
    """
    reg_list = []
    for chrom_name, chrom_len in chroms.items():
        regions = chr2reg(chrom_name, chrom_len, bin_size * 1000)
        if regions is None:
            print("Error: cannot split chromsome %s to bins!" % chrom_name)
            return None
        reg_list.extend(regions)
    return reg_list

def get_fixsize_reg_from_sam_header(chr_names, bin_size, sam_file):
    """
    @abstract         Create fixed-size bins for certain chroms acoording to chrom length in sam header
    @param chr_names  A list/tuple/array of chrom names to be splitted into bins [list/tuple/array]
    @param bin_size   Size of bin in kb [int]
    @param sam_file   Path to sam/bam/cram file [str]
    @return           A list of Region objects if success, None otherwise [list]
    """
    sam_fp = pysam.AlignmentFile(sam_file, "r")   # file format auto-detected.
    chroms = {}
    for chrom_name in chr_names:
        if chrom_name not in sam_fp.references:
            chrom_name = chrom_name[3:] if chrom_name.startswith("chr") else "chr" + chrom_name
            if chrom_name not in sam_fp.references:
                return None
        chrom_len = sam_fp.get_reference_length(chrom_name)
        chroms[chrom_name] = chrom_len
    sam_fp.close()
    return get_fixsize_reg_from_input_len(chroms, bin_size)

def get_fixsize_regions(bin_size, hg_ver):
    """
    @abstract        Create fixed-size bins for whole genome.
    @param bin_size  Size of bin in kb [int]
    @param hg_ver    Version of human genome: 19 or 38 [int]
    @return          A list of Region objects if success, None otherwise [list]
    """
    try:
        hg_ver = int(hg_ver)
    except:
        return None
    chr_names = [str(i) for i in range(1, 23)] + ["X", "Y"]
    chroms = {}
    if hg_ver == 19:
        chroms = {k:v for k, v in zip(chr_names, CHROM_LEN_HG19)}
    elif hg_ver == 38:
        chroms = {k:v for k, v in zip(chr_names, CHROM_LEN_HG38)}
    else:
        return None
    return get_fixsize_reg_from_input_len(chroms, bin_size)

def output_regions(reg_list, fname, ftype):
    """
    @abstract        Output Region objects to file
    @param reg_list  A list of Region objects to be outputed [list]
    @param fname     Path to the output file, use None for stdout [str]
    @param ftype     File type, one of bed|tsv [str]
    @return          0 if success, negative integer if error [int]
    """
    fp = None
    if fname is None: fp = sys.stdout
    elif fname: fp = open(fname, "w")
    else: return -1

    if ftype == "bed":
        for reg in reg_list:
            rec = "\t".join([reg.chrom, str(reg.start - 1), str(reg.stop), reg.id])
            fp.write(rec + "\n")
    elif ftype == "tsv":
        for reg in reg_list:
            rec = "\t".join([reg.chrom, str(reg.start), str(reg.stop)])
            fp.write(rec + "\n")
    else:
        return -3
    
    if fname is not None:
        fp.close()
    return 0

def reg2bed(reg_list, bed_file):
    """
    @abstract        Output Region objects to bed file
    @param reg_list  A list of Region objects to be outputed [list]
    @param bed_file  Path to the output bed file, use None for stdout [str]
    @return          0 if success, negative integer if error [int]
    """
    return output_regions(reg_list, bed_file, "bed")

def reg2tsv(reg_list, tsv_file):
    """
    @abstract        Output Region objects to tsv file
    @param reg_list  A list of Region objects to be outputed [list]
    @param tsv_file  Path to the output tsv file, use None for stdout [str]
    @return          0 if success, negative integer if error [int]
    """
    return output_regions(reg_list, tsv_file, "tsv")

