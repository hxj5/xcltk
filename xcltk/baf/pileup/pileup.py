# Aim: count spliced or un-spliced transcripts / reads in single cell data
#      to prepare inputs for RNA velocity modelling.
# Author: Xianjie Huang

# Note: now it only supports 10x data.

# TODO:
#   1. Add --maxDEPTH
#   2. Use pileup method besides fetch method

import getopt
import multiprocessing
import os
import pickle
import sys
import time

from .app import APP, VERSION
from .config import Config, \
                    CFG_DEBUG, \
                    CFG_CELL_TAG, CFG_UMI_TAG, CFG_UMI_TAG_BC, \
                    CFG_NPROC,     \
                    CFG_MIN_COUNT,  \
                    CFG_MIN_OVERLAP, CFG_MIN_UREAD, CFG_MIN_UFRAC, \
                    CFG_MIN_JUMI, CFG_MIN_JCELL,  \
                    CFG_DENOVO_JC_LEN, \
                    CFG_INCL_FLAG, CFG_EXCL_FLAG_UMI, CFG_EXCL_FLAG_XUMI, \
                    CFG_MIN_LEN, CFG_MIN_MAPQ
from .core import sp_count
from .gff import gff_load
from .thread import ThreadData
from .utils import merge_mtx, merge_tsv, rewrite_mtx
from .zfile import zopen, ZF_F_GZIP, ZF_F_PLAIN

def prepare_config(conf):
    """Prepare configures for downstream analysis
    @param conf  A Config object.
    @return      0 if success, -1 otherwise.
    @note        This function should be call after cmdline is parsed.
    """
    func = "prepare_config"
    conf.sam_fn_list = []
    #<DEV>#
    if conf.sam_fn:
        sam_lst0 = conf.sam_fn.split(",")
        for fn in sam_lst0:
            if os.path.isfile(fn):
                conf.sam_fn_list.append(fn)
            else:
                sys.stderr.write("[E::%s] failed to open sam file '%s'.\n" % (func, fn))
                return(-1)
    #elif conf.sam_list_fn:
    #    if os.path.isfile(conf.sam_list_fn):
    #        with zopen(conf.sam_list_fn, "rt") as fp:
    #            conf.sam_list = [pysam.AlignmentFile(x.strip(), "r") for x in fp]
    #    else:
    #        sys.stderr.write("[E::%s] failed to open sam list file '%s'.\n" % (func, conf.sam_list_fn))
    #        return(-1)
    else:
        sys.stderr.write("[E::%s] sam file(s) needed!\n" % (func,))
        return(-1)
    #</DEV>#

    if conf.barcode_fn:
        conf.sid_list = None
        if os.path.isfile(conf.barcode_fn):
            with zopen(conf.barcode_fn, "rt") as fp:
                conf.barcodes = sorted([x.strip() for x in fp])
            if len(set(conf.barcodes)) != len(conf.barcodes):
                sys.stderr.write("[E::%s] duplicate barcodes!\n" % (func, ))
                return(-1)
        else:
            sys.stderr.write("[E::%s] failed to open barcode file '%s'.\n" % (func, conf.barcode_fn))
            return(-1)
    else:
        #<DEV>#        
        conf.barcodes = None
        sys.stderr.write("[E::%s] barcode file needed!\n" % (func,))
        return(-1)
        #if conf.sid:
        #    conf.sid_list = [x.strip() for x in conf.sid.split(",")]
        #elif conf.sid_fn:
        #    if os.path.isfile(conf.sid_fn):
        #        with zopen(conf.sid_fn, "rt") as fp:
        #            conf.sid_list = [x.strip() for x in fp]
        #    else:
        #        sys.stderr.write("[E::%s] failed to open sample list file '%s'.\n" % (func, conf.sid_fn))
        #        return(-1)
        #else:
        #    conf.sid_list = ["Sample%d" % x for x in range(len(conf.sam_list))]
        #</DEV>#

    if not conf.out_dir:
        sys.stderr.write("[E::%s] out dir needed!\n" % func)
        return(-1)
    if not os.path.isdir(conf.out_dir):
        os.mkdir(conf.out_dir)
    conf.out_junction_fn = os.path.join(conf.out_dir, conf.out_prefix + "junction.tsv")
    conf.out_region_fn = os.path.join(conf.out_dir, conf.out_prefix + "region.tsv")
    conf.out_sample_fn = os.path.join(conf.out_dir, conf.out_prefix + "samples.tsv")
    conf.out_spliced_fn = os.path.join(conf.out_dir, conf.out_prefix + "spliced.mtx")
    conf.out_unspliced_fn = os.path.join(conf.out_dir, conf.out_prefix + "unspliced.mtx")
    conf.out_ambiguous_fn = os.path.join(conf.out_dir, conf.out_prefix + "ambiguous.mtx")

    if conf.region_fn:
        if os.path.isfile(conf.region_fn): 
            ret, lst = gff_load(conf.region_fn, gene_tag = "gene",
                tran_tag = "transcript,mRNA", exon_tag = "exon", verbose = False)
            if ret < 0:
                sys.stderr.write("[E::%s] failed to load region file.\n" % func)
                return(-1)
            conf.reg_list = lst.get_genes()
            if conf.barcodes is not None:
                sys.stdout.write("[I::%s] count %d regions in %d single cells.\n" % (func, 
                    len(conf.reg_list), len(conf.barcodes)))
            else:
                sys.stdout.write("[I::%s] count %d regions in %d bam files.\n" % (func,
                    len(conf.reg_list), len(conf.sid_list)))
        else:
            sys.stderr.write("[E::%s] failed to open region file '%s'.\n" % (func, conf.region_fn))
            return(-1)
    else:
        sys.stderr.write("[E::%s] region file needed!\n" % (func,))
        return(-1)

    if conf.cell_tag and conf.cell_tag.upper() == "NONE":
        conf.cell_tag = None
    if conf.cell_tag and conf.barcodes:
        if conf.sid_fn or conf.sid_list:
            sys.stderr.write("[E::%s] should not specify barcodes and sample IDs at the same time.\n" % (func, ))
            return(-1)
        conf.sid_list = None      
    elif (not conf.cell_tag) ^ (not conf.barcodes):
        sys.stderr.write("[E::%s] should not specify cell_tag or barcodes alone.\n" % (func, ))
        return(-1)
    else:
        #<DEV>#
        sys.stderr.write("[E::%s] should specify cell_tag and barcodes.\n" % (func, ))
        return(-1)
        #conf.cell_tag = conf.barcodes = None
        #if not conf.sid_fn or not conf.sid_list:
        #    sys.stderr.write("[E::%s] should specify either barcodes or sample IDs.\n" % (func, ))
        #    return(-1)
        #if len(conf.sid_list) != len(conf.sam_list):
        #    sys.stderr.write("[E::%s] len(sam_list) != len(sid_list).\n" % (func, ))
        #    return(-1)
        #</DEV>#         

    if conf.umi_tag:
        if conf.umi_tag.upper() == "AUTO":
            if conf.barcodes is None:
                conf.umi_tag = None
            else:
                conf.umi_tag = CFG_UMI_TAG_BC
        elif conf.umi_tag.upper() == "NONE":
            conf.umi_tag = None
    #<DEV>#
    if not conf.umi_tag:
        sys.stderr.write("[E::%s] umi tag needed!\n" % (func, ))
        return(-1)
    #</DEV>#

    if conf.barcodes:
        with open(conf.out_sample_fn, "w") as fp:
            fp.write("".join([b + "\n" for b in conf.barcodes]))

    return(0)

def usage(fp = sys.stderr):
    s =  "\n" 
    s += "Version: %s\n" % VERSION
    s += "Usage:   %s <options>\n" % APP  
    s += "\n" 
    s += "Options:\n"
    s += "  -s, --sam STR          Indexed sam/bam/cram file.\n"
    #<DEV/>#s += "  -s, --sam STR          Indexed sam/bam/cram file(s), comma separated multiple samples.\n"
    #<DEV/>#s += "  -S, --samList FILE     A list file containing bam files, each per line.\n"
    s += "  -O, --outdir DIR       Output directory for sparse matrices.\n"
    s += "  -R, --region FILE      A GTF/GFF3 file listing target regions.\n"
    s += "  -b, --barcode FILE     A plain file listing all effective cell barcode.\n"
    #<DEV/>#s += "  -i, --sampleIDs STR    Comma separated sample ids.\n"
    #<DEV/>#s += "  -I, --sampleLIST FILE  A list file containing sample IDs, each per line.\n"
    s += "  -V, --version          Print software version and exit.\n"
    s += "  -h, --help             Print this message and exit.\n"
    s += "  -D, --debug INT        Used by developer for debugging [%d]\n" % CFG_DEBUG
    s += "\n"
    s += "Optional arguments:\n"
    s += "  -p, --nproc INT        Number of processes [%d]\n" % CFG_NPROC
    #<DEV/>#s += "  --cellTAG STR          Tag for cell barcodes, turn off with None [%s]\n" % CFG_CELL_TAG
    s += "  --cellTAG STR          Tag for cell barcodes [%s]\n" % CFG_CELL_TAG
    #<DEV/>#s += "  --UMItag STR           Tag for UMI: TAG, Auto, None. For Auto mode, use %s if barcodes is inputted,\n"
    #<DEV/>#s += "                         otherwise use None. None mode means no UMI but read counts [%s]\n" % (CFG_UMI_TAG_BC, CFG_UMI_TAG)
    s += "  --UMItag STR           Tag for UMI [%s]\n" % CFG_UMI_TAG
    #s += "  --minCOUNT INT         Mininum aggragated count [%d]\n" % CFG_MIN_COUNT
    #<DEV/>#s += "  --pairEND              If use, analysis in pair-end mode.\n"
    s += "\n"
    s += "Junction parameters:\n"
    s += "  --minOVERLAP INT       Minimum overlapping bases with one region [%d]\n" % CFG_MIN_OVERLAP
    s += "  --minUREAD INT         Minimum number of supporting reads within one UMI [%d]\n" % CFG_MIN_UREAD
    s += "  --minUFRAC FLOAT       Minimum fraction of supporting reads within one UMI [%f]\n" % CFG_MIN_UFRAC
    s += "  --minJUMI INT          Minimum number of supporting UMIs for one junction [%d]\n" % CFG_MIN_JUMI
    s += "  --minJCELL INT         Minimum number of supporting cells for one junction [%d]\n" % CFG_MIN_JCELL
    #<DEV/>#s += "  --inclDENOVO           If set, include denovo junctions for analysis.\n"
    #<DEV/>#s += "  --denovoJCLen INT|FLOAT  Length (INT) or quantile (FLOAT) for selecting denovo junctions [%f]\n" % CFG_DENOVO_JC_LEN
    #<DEV/># add option to specify the maximum fraction of UMI that cross one denovo junction. (avoid calling germline large indel).
    s += "\n"
    s += "Read filtering:\n"
    s += "  --inclFLAG INT    Required flags: skip reads with all mask bits unset [%d]\n" % CFG_INCL_FLAG
    s += "  --exclFLAG INT    Filter flags: skip reads with any mask bits set [%d\n" % CFG_EXCL_FLAG_UMI
    s += "                    (when use UMI) or %d (otherwise)]\n" % CFG_EXCL_FLAG_XUMI
    s += "  --minLEN INT      Minimum mapped length for read filtering [%d]\n" % CFG_MIN_LEN
    s += "  --minMAPQ INT     Minimum MAPQ for read filtering [%d]\n" % CFG_MIN_MAPQ
    s += "  --countORPHAN     If use, do not skip anomalous read pairs.\n"
    s += "\n"

    fp.write(s)

def show_progress(rv = None):
    return(rv)

def main(argv):
    """Core part
    @param argv   A list of cmdline parameters [list]
    @return       0 if success, -1 otherwise [int]
    """
    func = "main"
    ret = -1

    start_time = time.time()

    conf = Config()
    opts, args = getopt.getopt(argv[1:], "-s:-S:-O:-R:-b:-i:-I:-V-h-D:-p:", [
                     "sam=", "samList=", 
                     "outdir=", 
                     "region=", "barcode=", 
                     "sampleIDs=", "sampleLIST=", 
                     "version", "help", "debug=",
                     "nproc=", 
                     "cellTAG=", "UMItag=", 
                     "minCOUNT=", "pairEND",
                     "minOVERLAP=", "minUREAD=", "minUFRAC=",
                     "minJUMI=", "minJCELL=",
                     "inclDENOVO", "denovoJCLen=",
                     "inclFLAG=", "exclFLAG=", "minLEN=", "minMAPQ=", "countORPHAN"
                ])

    for op, val in opts:
        if len(op) > 2:
            op = op.lower()
        if op in   ("-s", "--sam"): conf.sam_fn = val
        elif op in ("-S", "--samlist"): conf.sam_list_fn = val
        elif op in ("-O", "--outdir"): conf.out_dir = val
        elif op in ("-R", "--region"): conf.region_fn = val
        elif op in ("-b", "--barcode"): conf.barcode_fn = val
        elif op in ("-i", "--sampleids"): conf.sid = val
        elif op in ("-I", "--samplelist"): conf.sid_fn = val
        elif op in ("-V", "--version"): sys.stderr.write(VERSION + "\n"); sys.exit(1)
        elif op in ("-h", "--help"): usage(sys.stderr); sys.exit(1)
        elif op in ("-D", "--debug"): conf.debug = int(val)

        elif op in ("-p", "--proc"): conf.nproc = int(val)
        elif op in ("--celltag"): conf.cell_tag = val
        elif op in ("--umitag"): conf.umi_tag = val
        elif op in ("--mincount"): conf.min_count = int(val)
        elif op in ("--pairend"): conf.pair_end = True

        elif op in ("--minoverlap"): conf.min_overlap = int(val)
        elif op in ("--minuread"): conf.min_uread = int(val)
        elif op in ("--minufrac"): conf.min_ufrac = float(val)
        elif op in ("--minjumi"): conf.min_jumi = int(val)
        elif op in ("--minjcell"): conf.min_jcell = int(val)
        elif op in ("--incldenovo"): conf.incl_denovo = True
        elif op in ("--denovojclen"): conf.denovo_jc_len = float(val)

        elif op in ("--inclflag"): conf.incl_flag = int(val)
        elif op in ("--exclflag"): conf.excl_flag = int(val)
        elif op in ("--minlen"): conf.min_len = int(val)
        elif op in ("--minmapq"): conf.min_mapq = float(val)
        elif op in ("--countorphan"): conf.no_orphan = False

        else:
            sys.stderr.write("[E::%s] invalid option: '%s'.\n" % (func, op))
            return(-1)

    try:
        time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(start_time))
        sys.stdout.write("[I::%s] start time: %s.\n" % (func, time_str))

        cmdline = " ".join(argv)
        sys.stdout.write("[I::%s] CMD: %s\n" % (func, cmdline))
        sys.stdout.write("[I::%s] VERSION: %s %s\n" % (func, APP, VERSION))

        if prepare_config(conf) < 0:
            raise ValueError("[%s] errcode %d" % (func, -2))
        sys.stderr.write("[D::%s] program configuration:\n" % func)
        conf.show(fp = sys.stderr, prefix = "\t")

        # For multiprocessing, we can split the gene set in two ways to save memory:
        # 1. implement a generator (e.g., with 'yield') to push one gene to Process pool every time;
        #    and merge & sort all genes in the end.
        # 2. read all genes from gtf/gff file -> split gene set -> save each sub-gene-list to disk
        #    (e.g., in pickle format) -> read each sub-gene-list in one sub-Process.
        # It seems strategy 2 is a better choice (easier to implement).
        
        m_reg = len(conf.reg_list)
        m_thread = conf.nproc if m_reg >= conf.nproc else m_reg

        reg_fn_list = []
        if m_thread > 1:
            n_reg = m_reg // m_thread
            r_reg = m_reg - n_reg * m_thread
            k_reg = 0
            i_thread = 0
            while k_reg <= m_reg - 1:
                t_reg = n_reg + 1 if i_thread < r_reg else n_reg
                reg_fn = conf.out_prefix + "region.pickle." + str(i_thread)
                reg_fn = os.path.join(conf.out_dir, reg_fn)
                reg_fn_list.append(reg_fn)
                with open(reg_fn, "wb") as fp:
                    pickle.dump(conf.reg_list[k_reg:(k_reg + t_reg)], fp)
                k_reg += t_reg
                i_thread += 1
            for gen in conf.reg_list:  # save memory
                del gen
            conf.reg_list.clear()
            conf.reg_list = None

        thdata_list = []
        ret_sp = -1
        if m_thread <= 1:
            thdata = ThreadData(
                idx = 0, conf = conf,
                reg_obj = conf.reg_list, is_reg_pickle = False,
                out_junction_fn = conf.out_junction_fn + ".0",
                out_region_fn = conf.out_region_fn + ".0",
                out_spliced_fn = conf.out_spliced_fn + ".0",
                out_unspliced_fn = conf.out_unspliced_fn + ".0",
                out_ambiguous_fn = conf.out_ambiguous_fn + ".0",
                out_fn = None
            )
            thdata_list.append(thdata)
            if conf.debug > 0:
                sys.stderr.write("[D::%s] data of thread-%d before sp_count:\n" % (func, 0))
                thdata.show(fp = sys.stderr, prefix = "\t")
            ret_sp, thdata = sp_count(thdata)
        else:
            pool = multiprocessing.Pool(processes = m_thread)
            mp_result = []
            for i in range(m_thread):
                thdata = ThreadData(
                    idx = i, conf = conf,
                    reg_obj = reg_fn_list[i], is_reg_pickle = True,
                    out_junction_fn = conf.out_junction_fn + "." + str(i),
                    out_region_fn = conf.out_region_fn + "." + str(i),
                    out_spliced_fn = conf.out_spliced_fn + "." + str(i),
                    out_unspliced_fn = conf.out_unspliced_fn + "." + str(i),
                    out_ambiguous_fn = conf.out_ambiguous_fn + "." + str(i),
                    out_fn = None   
                )
                thdata_list.append(thdata)
                if conf.debug > 0:
                    sys.stderr.write("[D::%s] data of thread-%d before sp_count:\n" % (func, i))
                    thdata.show(fp = sys.stderr, prefix = "\t")
                mp_result.append(pool.apply_async(
                        func = sp_count, 
                        args = (thdata, ), 
                        callback = show_progress))   # TODO: error_callback?
            pool.close()
            pool.join()
            mp_result = [res.get() for res in mp_result]
            retcode_list = [item[0] for item in mp_result]
            thdata_list = [item[1] for item in mp_result]
            if conf.debug > 0:
                sys.stderr.write("[D::%s] returned values of multi-processing:\n" % func)
                sys.stderr.write("\t%s\n" % str(retcode_list))

        # check running status of each sub-process
        for thdata in thdata_list:         
            if conf.debug > 0:
                sys.stderr.write("[D::%s] data of thread-%d after sp_count:\n" % (func, thdata.idx))
                thdata.show(fp = sys.stderr, prefix = "\t")
            if thdata.ret < 0:
                raise ValueError("[%s] errcode %d" % (func, -3))

        if m_thread <= 1:
            if ret_sp < 0:
                raise ValueError("[%s] errcode %d" % (func, -5))

            thdata = thdata_list[0]
            os.rename(thdata.out_junction_fn, conf.out_junction_fn)
            os.rename(thdata.out_region_fn, conf.out_region_fn)

            if rewrite_mtx(thdata.out_spliced_fn, ZF_F_GZIP, 
                           conf.out_spliced_fn, "wb", ZF_F_PLAIN, 
                           m_reg, len(conf.barcodes), thdata.nr_sp,
                           remove = True) < 0:
                raise ValueError("[%s] errcode %d" % (func, -7))

            if rewrite_mtx(thdata.out_unspliced_fn, ZF_F_GZIP, 
                           conf.out_unspliced_fn, "wb", ZF_F_PLAIN, 
                           m_reg, len(conf.barcodes), thdata.nr_us,
                           remove = True) < 0:
                raise ValueError("[%s] errcode %d" % (func, -9)) 

            if rewrite_mtx(thdata.out_ambiguous_fn, ZF_F_GZIP, 
                           conf.out_ambiguous_fn, "wb", ZF_F_PLAIN, 
                           m_reg, len(conf.barcodes), thdata.nr_am,
                           remove = True) < 0:
                raise ValueError("[%s] errcode %d" % (func, -11))
        else:
            if merge_tsv([td.out_junction_fn for td in thdata_list], ZF_F_GZIP,
                         conf.out_junction_fn, "wb", ZF_F_PLAIN,
                         remove = True) < 0:
                raise ValueError("[%s] errcode %d" % (func, -13))

            if merge_tsv([td.out_region_fn for td in thdata_list], ZF_F_GZIP, 
                         conf.out_region_fn, "wb", ZF_F_PLAIN, 
                         remove = True) < 0:
                raise ValueError("[%s] errcode %d" % (func, -15))

            nr_reg_list = [td.nr_reg for td in thdata_list]

            if merge_mtx([td.out_spliced_fn for td in thdata_list], ZF_F_GZIP, 
                         conf.out_spliced_fn, "w", ZF_F_PLAIN,
                         nr_reg_list, len(conf.barcodes), sum([td.nr_sp for td in thdata_list]),
                         remove = True) < 0:
                raise ValueError("[%s] errcode %d" % (func, -17))

            if merge_mtx([td.out_unspliced_fn for td in thdata_list], ZF_F_GZIP, 
                         conf.out_unspliced_fn, "w", ZF_F_PLAIN,
                         nr_reg_list, len(conf.barcodes), sum([td.nr_us for td in thdata_list]),
                         remove = True) < 0:
                raise ValueError("[%s] errcode %d" % (func, -19))

            if merge_mtx([td.out_ambiguous_fn for td in thdata_list], ZF_F_GZIP, 
                         conf.out_ambiguous_fn, "w", ZF_F_PLAIN,
                         nr_reg_list, len(conf.barcodes), sum([td.nr_am for td in thdata_list]),
                         remove = True) < 0:
                raise ValueError("[%s] errcode %d" % (func, -21))

    except ValueError as e:
        sys.stderr.write("[E::%s] '%s'\n" % (func, str(e)))
        sys.stdout.write("[E::%s] Running program failed.\n" % func)
        sys.stdout.write("[E::%s] Quiting ...\n" % func)
        ret = -1

    else:
        sys.stdout.write("[I::%s] All Done!\n" % func)
        ret = 0

    finally:
        sys.stdout.write("[I::%s] CMD: %s\n" % (func, cmdline))
        sys.stdout.write("[I::%s] VERSION: %s %s\n" % (func, APP, VERSION))

        end_time = time.time()
        time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(end_time))
        sys.stdout.write("[I::%s] end time: %s\n" % (func, time_str))
        sys.stdout.write("[I::%s] time spent: %.2fs\n" % (func, end_time - start_time))

    return(ret)

def run():
    if len(sys.argv) < 2:
        usage(sys.stderr)
        sys.exit(1)

    main(sys.argv)

if __name__ == "__main__":
    run()
