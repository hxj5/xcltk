# main.py - basic feature counting.


import getopt
import multiprocessing
import os
import pickle
import sys
import time

from logging import debug, error, info
from logging import warning as warn

from .config import Config
from .core import fc_features
from .thread import ThreadData
from .utils import load_region_from_txt, merge_mtx, merge_tsv

from ...config import APP, VERSION
from ...utils.xlog import init_logging
from ...utils.zfile import zopen, ZF_F_GZIP, ZF_F_PLAIN

COMMAND = "basefc"


def usage(fp = sys.stdout, conf = None):
    s =  "\n"
    s += "Version: %s\n" % VERSION
    s += "Usage:   %s %s <options>\n" % (APP, COMMAND)
    s += "\n" 
    s += "Options:\n"
    s += "  -s, --sam FILE         Comma separated indexed sam/bam/cram file.\n"
    s += "  -S, --samList FILE     A list file containing bam files, each per line.\n"
    s += "  -b, --barcode FILE     A plain file listing all effective cell barcode.\n"
    s += "  -R, --region FILE      A TSV file listing target regions. The first 4 columns shoud be:\n"
    s += "                         chrom, start, end (both 1-based and inclusive), name.\n"
    s += "  -i, --sampleList FILE  A list file containing sample IDs, each per line.\n"
    s += "  -I, --sampleIDs STR    Comma separated sample IDs.\n"
    s += "  -O, --outdir DIR       Output directory for sparse matrices.\n"
    s += "  -h, --help             Print this message and exit.\n"
    s += "\n"
    s += "Optional arguments:\n"
    s += "  -p, --ncores INT       Number of processes [%d]\n" % conf.NPROC
    s += "      --cellTAG STR      Tag for cell barcodes, set to None when using sample IDs [%s]\n" % conf.CELL_TAG
    s += "      --UMItag STR       Tag for UMI, set to None when reads only [%s]\n" % conf.UMI_TAG
    s += "  -D, --debug INT        Used by developer for debugging [%d]\n" % conf.DEBUG
    s += "\n"
    s += "Read filtering:\n"
    s += "  --inclFLAG INT          Required flags: skip reads with all mask bits unset [%d]\n" % conf.INCL_FLAG
    s += "  --exclFLAG INT          Filter flags: skip reads with any mask bits set [%d\n" % conf.EXCL_FLAG_UMI
    s += "                          (when use UMI) or %d (otherwise)]\n" % conf.EXCL_FLAG_XUMI
    s += "  --minLEN INT            Minimum mapped length for read filtering [%d]\n" % conf.MIN_LEN
    s += "  --minMAPQ INT           Minimum MAPQ for read filtering [%d]\n" % conf.MIN_MAPQ
    s += "  --minINCLUDE FLOAT|INT  Minimum fraction or length of included part within specific feature [%f]\n" % conf.MIN_INCLUDE
    s += "  --countORPHAN           If use, do not skip anomalous read pairs.\n"
    s += "\n"

    fp.write(s)


def fc_main(argv, conf = None):
    """Command-Line interface.

    Parameters
    ----------
    argv : list
        A list of cmdline parameters.
    conf : fc::Config object
        Configuration object.
    
    Returns
    -------
    int
        0 if success, -1 otherwise [int]
    """
    if conf is None:
        conf = Config()

    if len(argv) <= 2:
        usage(sys.stdout, conf.defaults)
        sys.exit(0)

    conf.argv = argv.copy()
    init_logging(stream = sys.stderr)

    opts, args = getopt.getopt(
        args = argv[2:], 
        shortopts = "-s:-S:-b:-R:-i:-I:-O:-h-p:-D:",
        longopts = [
            "sam=", "samList=", "barcode=",
            "region=",
            "sampleList=", "sampleIDs=",
            "outdir=",
            "help",

            "ncores=",
            "cellTAG=", "UMItag=",
            "debug=",

            "inclFLAG=", "exclFLAG=", 
            "minLEN=", "minMAPQ=", 
            "minINCLUDE=",
            "countORPHAN"
        ])

    for op, val in opts:
        if len(op) > 2:
            op = op.lower()
        if op in   ("-s", "--sam"): conf.sam_fn = val
        elif op in ("-S", "--samlist"): conf.sam_list_fn = val
        elif op in ("-b", "--barcode"): conf.barcode_fn = val
        elif op in ("-R", "--region"): conf.region_fn = val
        elif op in ("-i", "--samplelist"): conf.sample_id_fn = val
        elif op in ("-I", "--sampleids"): conf.sample_id_str = val
        elif op in ("-O", "--outdir"): conf.out_dir = val
        elif op in ("-h", "--help"): usage(sys.stdout, conf.defaults); sys.exit(0)

        elif op in ("-p", "--ncores"): conf.nproc = int(val)
        elif op in (      "--celltag"): conf.cell_tag = val
        elif op in (      "--umitag"): conf.umi_tag = val
        elif op in ("-D", "--debug"): conf.debug = int(val)

        elif op in ("--inclflag"): conf.incl_flag = int(val)
        elif op in ("--exclflag"): conf.excl_flag = int(val)
        elif op in ("--minlen"): conf.min_len = int(val)
        elif op in ("--minmapq"): conf.min_mapq = float(val)
        elif op in ("--mininclude"):
            if "." in val:
                conf.min_include = float(val)
            else:
                conf.min_include = int(val)
        elif op in ("--countorphan"): conf.no_orphan = False

        else:
            error("invalid option: '%s'." % op)
            return(-1)
        
    ret = fc_run(conf)
    return(ret)


def fc_wrapper(
    sam_fn, barcode_fn,
    region_fn,
    out_dir,
    sam_list_fn = None,
    sample_ids = None, sample_id_fn = None,
    debug_level = 0,
    ncores = 1,
    cell_tag = "CB", umi_tag = "UB",
    output_all_reg = True,
    min_mapq = 20, min_len = 30,
    min_include = 0.9,
    incl_flag = 0, excl_flag = None,
    no_orphan = True
):
    conf = Config()

    conf.sam_fn = sam_fn
    conf.sam_list_fn = sam_list_fn
    conf.barcode_fn = barcode_fn
    conf.region_fn = region_fn
    conf.sample_id_str = sample_ids
    conf.sample_id_fn = sample_id_fn
    conf.out_dir = out_dir
    conf.debug = debug_level

    conf.cell_tag = cell_tag
    conf.umi_tag = umi_tag
    conf.nproc = ncores
    conf.output_all_reg = output_all_reg

    conf.min_mapq = min_mapq
    conf.min_len = min_len
    conf.min_include = min_include
    conf.incl_flag = incl_flag
    if excl_flag is None:
        conf.excl_flag = -1
    conf.no_orphan = no_orphan

    ret = fc_run(conf)
    return(ret)


def fc_core(conf):
    if prepare_config(conf) < 0:
        raise ValueError("errcode -2")
    info("program configuration:")
    conf.show(fp = sys.stderr, prefix = "\t")

    # split region list and save to file
    m_reg = len(conf.reg_list)
    m_thread = conf.nproc if m_reg >= conf.nproc else m_reg

    reg_fn_list = []
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
    for reg in conf.reg_list:  # save memory
        del reg
    conf.reg_list.clear()
    conf.reg_list = None

    thdata_list = []
    pool = multiprocessing.Pool(processes = m_thread)
    mp_result = []
    for i in range(m_thread):
        thdata = ThreadData(
            idx = i, conf = conf,
            reg_obj = reg_fn_list[i], is_reg_pickle = True,
            out_region_fn = conf.out_region_fn + "." + str(i),
            out_mtx_fn = conf.out_mtx_fn + "." + str(i),
            out_fn = None
        )
        thdata_list.append(thdata)
        if conf.debug > 0:
            debug("data of thread-%d before fc_features:" % i)
            thdata.show(fp = sys.stderr, prefix = "\t")
        mp_result.append(pool.apply_async(
            func = fc_features, 
            args = (thdata, ), 
            callback = show_progress))   # TODO: error_callback?
    pool.close()
    pool.join()
    mp_result = [res.get() for res in mp_result]
    retcode_list = [item[0] for item in mp_result]
    thdata_list = [item[1] for item in mp_result]
    if conf.debug > 0:
        debug("returned values of multi-processing:")
        debug("\t%s" % str(retcode_list))

    # check running status of each sub-process
    for thdata in thdata_list:         
        if conf.debug > 0:
            debug("data of thread-%d after fc_features:" %  thdata.idx)
            thdata.show(fp = sys.stderr, prefix = "\t")
        if thdata.ret < 0:
            raise ValueError("errcode -3")

    # merge results
    if merge_tsv(
        [td.out_region_fn for td in thdata_list], ZF_F_GZIP, 
        conf.out_region_fn, "wb", ZF_F_PLAIN, 
        remove = True
    ) < 0:
        raise ValueError("errcode -15")

    nr_reg_list = [td.nr_reg for td in thdata_list]

    if merge_mtx(
        [td.out_mtx_fn for td in thdata_list], ZF_F_GZIP,
        conf.out_mtx_fn, "w", ZF_F_PLAIN,
        nr_reg_list, len(conf.samples),
        sum([td.nr_mtx for td in thdata_list]),
        remove = True
    ) < 0:
        raise ValueError("errcode -17")
    

def fc_run(conf):
    ret = -1
    cmdline = None

    start_time = time.time()
    time_str = time.strftime(
        "%Y-%m-%d %H:%M:%S", time.localtime(start_time))
    info("start time: %s." % time_str)

    if conf.argv is not None:
        cmdline = " ".join(conf.argv)
        info("CMD: %s" % cmdline)

    try:
        ret = fc_core(conf)
    except ValueError as e:
        error(str(e))
        error("Running program failed.")
        error("Quiting ...")
        ret = -1
    else:
        info("All Done!")
        ret = 0
    finally:
        if conf.argv is not None:
            info("CMD: %s" % cmdline)

        end_time = time.time()
        time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(end_time))
        info("end time: %s" % time_str)
        info("time spent: %.2fs" % (end_time - start_time, ))

    return(ret)
    

def prepare_config(conf):
    """Prepare configures for downstream analysis

    Parameters
    ----------
    conf :  Config object
        Configuration info.

    Returns
    -------
    int
        0 if success, -1 otherwise.

    Notes
    -----
    This function should be called after cmdline is parsed.
    """
    if conf.sam_fn:
        if conf.sam_list_fn:
            error("should not specify 'sam_fn' and 'sam_list_fn' together.")
            return(-1)
        conf.sam_fn_list = conf.sam_fn.split(",")
    else:
        if not conf.sam_list_fn:
            error("one of 'sam_fn' and 'sam_list_fn' should be specified.")
            return(-1)
        with open(conf.sam_list_fn, "r") as fp:
            conf.sam_fn_list = [x.rstrip() for x in fp.readlines()]
    
    for fn in conf.sam_fn_list:
        if not os.path.isfile(fn):
            error("sam file '%s' does not exist." % fn)
            return(-1)

    if conf.barcode_fn:
        conf.sample_ids = None
        if conf.sample_id_str or conf.sample_id_fn:
            error("should not specify barcodes and sample IDs together.")
            return(-1)
        if os.path.isfile(conf.barcode_fn):
            with zopen(conf.barcode_fn, "rt") as fp:
                conf.barcodes = sorted([x.strip() for x in fp])   # UPDATE!! use numpy or pandas to load
            if len(set(conf.barcodes)) != len(conf.barcodes):
                error("duplicate barcodes!")
                return(-1)
        else:
            error("barcode file '%s' does not exist." % conf.barcode_fn)
            return(-1)
    else:
        conf.barcodes = None
        if conf.sample_id_str and conf.sample_id_fn:
            error("should not specify 'sample_id_str' and 'sample_fn' together.")
            return(-1)
        elif conf.sample_id_str:
            conf.sample_ids = conf.sample_id_str.split(",")
        elif conf.sample_id_fn:
            with zopen(conf.sample_id_fn, "rt") as fp:
                conf.sample_ids = [x.strip() for x in fp]
        else:
            warn("use default sample IDs ...")
            conf.sample_ids = ["Sample%d" % i for i in \
                range(len(conf.sam_fn_list))]
        if len(conf.sample_ids) != len(conf.sam_fn_list):
            error("numbers of sam files and sample IDs are different.")
            return(-1)
        
    conf.samples = conf.barcodes if conf.barcodes else conf.sample_ids

    if not conf.out_dir:
        error("out dir needed!")
        return(-1)
    if not os.path.isdir(conf.out_dir):
        os.mkdir(conf.out_dir)
    conf.out_region_fn = os.path.join(
        conf.out_dir, conf.out_prefix + "features.tsv")
    conf.out_sample_fn = os.path.join(
        conf.out_dir, conf.out_prefix + "barcodes.tsv")
    conf.out_mtx_fn = os.path.join(
        conf.out_dir, conf.out_prefix + "matrix.mtx")

    if conf.region_fn:
        if os.path.isfile(conf.region_fn): 
            conf.reg_list = load_region_from_txt(
                conf.region_fn, verbose = True)
            if not conf.reg_list:
                error("failed to load region file.")
                return(-1)
            info("count %d regions in %d single cells." % (
                len(conf.reg_list), len(conf.samples)))
        else:
            error("region file '%s' does not exist." % conf.region_fn)
            return(-1)
    else:
        error("region file needed!")
        return(-1)

    if conf.cell_tag and conf.cell_tag.upper() == "NONE":
        conf.cell_tag = None
    if conf.cell_tag and conf.barcodes:
        pass       
    elif (not conf.cell_tag) ^ (not conf.barcodes):
        error("should not specify cell_tag or barcodes alone.")
        return(-1)
    else:
        pass    

    if conf.umi_tag:
        if conf.umi_tag.upper() == "AUTO":
            if conf.barcodes is None:
                conf.umi_tag = None
            else:
                conf.umi_tag = conf.defaults.UMI_TAG_BC
        elif conf.umi_tag.upper() == "NONE":
            conf.umi_tag = None
    else:
        pass

    with open(conf.out_sample_fn, "w") as fp:
        fp.write("".join([smp + "\n" for smp in conf.samples]))

    if conf.excl_flag < 0:
        if conf.use_umi():
            conf.excl_flag = conf.defaults.EXCL_FLAG_UMI
        else:
            conf.excl_flag = conf.defaults.EXCL_FLAG_XUMI

    return(0)


def show_progress(rv = None):
    return(rv)
