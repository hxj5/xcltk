# Count reads for features in each cell
# Author: Yuanhua Huang
# Date: 27-04-2020
# Modified by Xianjie Huang 

import os
import sys
import time
import numpy as np
import multiprocessing
from optparse import OptionParser, OptionGroup

from .config import APP
from ..config import VERSION
from ..utils.count import feature_count
from ..utils.region import load_regions

FID = None
PROCESSED = 0
TOTAL_REGION = 0
START_TIME = time.time()

COMMAND = "basefc"

def show_progress(RV=None):
    global PROCESSED, TOTAL_REGION, START_TIME, FID
    if RV is not None:
        FID.writelines(RV)
    
    PROCESSED += 1
    bar_len = 20
    run_time = time.time() - START_TIME
    percents = 100.0 * PROCESSED / TOTAL_REGION
    filled_len = int(round(bar_len * percents / 100))
    bar = '=' * filled_len + '-' * (bar_len - filled_len)
    
    sys.stderr.write('[%s] [%s] %.1f%% done in %.1f sec.\n' 
        % (COMMAND, bar, percents, run_time))
    sys.stderr.flush()
    return RV

def base_fc(argv):
    # import warnings
    # warnings.filterwarnings('error')

    # parse command line options
    if len(argv) < 3:
        print("Welcome to %s %s v%s!\n" % (APP, COMMAND, VERSION))
        print("use -h or --help for help on argument.")
        sys.exit(1)
        
    parser = OptionParser(usage = "Usage: %s %s [options]" % (APP, COMMAND))
    parser.add_option("--samFile", "-s", dest="sam_file", default=None,
        help=("Indexed & sorted sam/bam/cram file, possibly from CellRanger."))
    parser.add_option("--outDir", "-O", dest="out_dir", default=None,
        help=("Full path for output direcotry [default: $samFile_dir]"))
    parser.add_option("--barcodeFile", "-b", dest="barcode_file", default=None,
        help=("A plain file listing all effective cell barcodes."))
    parser.add_option("--regFile", "-r", dest="reg_file", default=None,
        help=("Path to region file."))
    parser.add_option("--regType", "-T", dest="reg_type", default=None,
        help=("Region type, one of bed|gff|tsv; " 
              "tsv has 3 tab-delimed columns: <chr> <start> <end> "
              "and both start and end 1-based and included."))
    
    group1 = OptionGroup(parser, "Optional arguments")
    group1.add_option("--nproc", "-p", type="int", dest="nproc", default=1,
        help="Number of subprocesses [default: %default]")
    group1.add_option("--cellTAG", dest="cell_tag", default="CB", 
        help="Tag for cell barcodes, turn off with None [default: %default]")
    group1.add_option("--UMItag", dest="UMI_tag", default="None", 
        help="Tag for UMI: UR (common for 10x RNA). None means no UMI but read "
        "counts [default: %default]")
    
    group2 = OptionGroup(parser, "Read filtering")
    group2.add_option("--minLEN", type="int", dest="min_LEN", default=30, 
        help="Minimum length mapped to the feature [default: %default]")
    group2.add_option("--minMAPQ", type="int", dest="min_MAPQ", default=20, 
        help="Minimum MAPQ for read filtering [default: %default]")
    group2.add_option("--maxFLAG", type="int", dest="max_FLAG", default=255, 
        help="Maximum FLAG for read filtering [default: %default]")
    
    parser.add_option_group(group1)
    parser.add_option_group(group2)

    # check options and args
    (options, args) = parser.parse_args(args = argv[2:])
    if options.sam_file is None or options.barcode_file is None:
        print("Error: need samFile and barcodeFile.")
        sys.exit(1)

    sam_file = options.sam_file
    barcodes = np.genfromtxt(options.barcode_file, dtype="str", delimiter="\t")
    barcodes = sorted(list(barcodes))
    cell_tag = options.cell_tag
    UMI_tag  = options.UMI_tag if options.UMI_tag.upper() != "NONE" else None
        
    if options.out_dir is None:
        out_dir = os.path.dirname(options.sam_file)
    elif os.path.abspath(options.out_dir) == "":
        out_dir = "./" + options.out_dir
    else:
        out_dir = options.out_dir
    if os.path.isdir(os.path.abspath(out_dir)) == False:
        print("Error: No such directory for file\n -- %s" %out_dir)
        sys.exit(1)   

    reg_file, reg_type = options.reg_file, options.reg_type 
    if not reg_file or not reg_type:
        print("Error: both region file and region type should be provided!")
        sys.exit(1)
    elif not os.path.isfile(reg_file):
        print("Error: region file not exist: %s" % reg_file)
        sys.exit(1)
    elif reg_type.lower() not in ("bed", "gff", "tsv"):
        print("Error: region type should be one of bed|gff|tsv.")
        sys.exit(1)
    else:
        reg_type = reg_type.lower()

    nproc = options.nproc
    min_LEN = options.min_LEN
    min_MAPQ = options.min_MAPQ
    max_FLAG = options.max_FLAG
    
    # load regions
    reg_list = load_regions(reg_file, reg_type)
    if not reg_list:
        print("Error: empty region file or failed to parse region file!")
        sys.exit(1)

    reg_ids = [reg.id for reg in reg_list]
    feature_name = "gene" if reg_type == "gff" else "region"
    if len(np.unique(reg_ids)) != len(reg_ids):
        reg_ids = ["%s%d" % (feature_name, x + 1) for x in range(len(reg_list))]

    # save cell barcodes and regions
    fid = open(out_dir + "/features.tsv", "w")
    for _reg_id in reg_ids:
        fid.writelines(_reg_id + "\n")
    fid.close()
    print("[%s] Count reads for %d features in %d cells with %d cores." 
          % (COMMAND, len(reg_ids), len(barcodes), nproc))

    fid = open(out_dir + "/barcodes.tsv", "w")
    for _cell_id in barcodes:
        fid.writelines(_cell_id + "\n")
    fid.close()

    global FID, TOTAL_REGION
    TOTAL_REGION = len(reg_ids)
    FID = open(out_dir + "/matrix.mtx", "w")

    if nproc > 1:
        result = []
        pool = multiprocessing.Pool(processes=nproc)
        for ii in range(len(reg_list)):
            result.append(pool.apply_async(feature_count, (sam_file, 
                barcodes, reg_list[ii], ii, cell_tag, UMI_tag, min_MAPQ, max_FLAG, min_LEN), 
                callback=show_progress))
        pool.close()
        pool.join()
    else:
        for ii in range(len(reg_list)):
            RV = feature_count(sam_file, barcodes, reg_list[ii], ii, cell_tag, 
                               UMI_tag, min_MAPQ, max_FLAG, min_LEN)
            show_progress(RV)

    FID.close()

    data = np.genfromtxt(out_dir + "/matrix.mtx")
    fid = open(out_dir + "/matrix.mtx", "w")
    fid.writelines(
        "%" + "%MatrixMarket matrix coordinate integer general\n" + "%\n")
    fid.writelines("%d\t%d\t%d\n" % (len(reg_list), len(barcodes), data.shape[0]))
    
    sort_idx = np.argsort(data[:, 0] * len(barcodes) + data[:, 1])
    for ii in sort_idx:
        fid.writelines("%d\t%d\t%d\n" %(data[ii, 0] + 1, data[ii, 1] + 1, data[ii, 2]))  # set cell_idx and reg_idx to 1-based
    fid.close()
    
    print("")
    

if __name__ == "__main__":
    base_fc(sys.argv)

