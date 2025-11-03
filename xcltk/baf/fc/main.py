# main.py - allele-specific feature counting.


import gc
import getopt
import multiprocessing
import numpy as np
import os
import pickle
import sys
import time

from logging import debug, error, info
from logging import warning as warn

from .config import Config
from .core import fc_features
from .phasing import reg_local_phasing
from .thread import ThreadData
from .utils import load_region_from_txt, load_snp_from_vcf, \
    load_snp_from_tsv, merge_mtx, merge_tsv

from ...config import APP, VERSION
from ...utils.base import assert_e
from ...utils.csp_io import load_data as csp_load_data
from ...utils.grange import format_chrom
from ...utils.xlog import init_logging
from ...utils.zfile import zopen, ZF_F_GZIP, ZF_F_PLAIN



def afc_wrapper(
    sam_fn, barcode_fn,
    region_fn, phased_snp_fn, 
    out_dir,
    sam_list_fn = None,
    sample_ids = None, sample_id_fn = None,
    debug_level = 0,
    ncores = 1,
    cellsnp_dir = None, ref_cell_fn = None,
    cell_tag = "CB", umi_tag = "UB",
    min_count = 1, min_maf = 0,
    output_all_reg = False, no_dup_hap = True,
    min_mapq = 20, min_len = 30,
    incl_flag = 0, excl_flag = None,
    no_orphan = True
):
    conf = Config()

    conf.sam_fn = sam_fn
    conf.sam_list_fn = sam_list_fn
    conf.barcode_fn = barcode_fn
    conf.region_fn = region_fn
    conf.snp_fn = phased_snp_fn
    conf.sample_id_str = sample_ids
    conf.sample_id_fn = sample_id_fn
    conf.out_dir = out_dir
    conf.debug = debug_level

    conf.cellsnp_dir = cellsnp_dir
    conf.ref_cell_fn = ref_cell_fn
    conf.cell_tag = cell_tag
    conf.umi_tag = umi_tag
    conf.nproc = ncores
    conf.min_count = min_count
    conf.min_maf = min_maf
    conf.output_all_reg = output_all_reg
    conf.no_dup_hap = no_dup_hap

    conf.min_mapq = min_mapq
    conf.min_len = min_len
    conf.incl_flag = incl_flag
    conf.excl_flag = -1 if excl_flag is None else excl_flag
    conf.no_orphan = no_orphan

    ret = afc_run(conf)
    return(ret)



def afc_core(conf):
    if prepare_config(conf) < 0:
        raise ValueError("errcode -2")
    info("program configuration:")
    conf.show(fp = sys.stderr, prefix = "\t")

    
    # extract SNPs for each region
    if conf.debug > 0:
        debug("extract SNPs for each region.")
    reg_list = []
    for reg in conf.reg_list:
        snp_list = conf.snp_set.fetch(reg.chrom, reg.start, reg.end)
        if snp_list and len(snp_list) > 0:
            reg.snp_list = sorted(snp_list, key = lambda s: s.pos)
            reg_list.append(reg)
        else:
            if conf.debug > 2:
                debug("region '%s': no SNP fetched." % reg.name)
    info("#regions: total=%d; with_snps=%d." % \
         (len(conf.reg_list), len(reg_list)))

    if not conf.output_all_reg:
        conf.reg_list = reg_list
        
        
    # do local phasing within each region.
    n_rlp = 0             # regions that do local phasing.
    n_rlp_failed = 0
    n_slp = 0             # SNPs that do local phasing.
    n_slp_flipped = 0
    if conf.use_local_phasing():
        adata = conf.snp_adata
        if conf.ref_cells is not None:
            adata = adata[~adata.obs['cell'].isin(conf.ref_cells), :]
        adata.var['chrom'] = adata.var['chrom'].map(format_chrom)
        for reg in conf.reg_list:
            if reg.end - reg.start < conf.rlp_min_len:
                continue
            if reg.snp_list is None or len(reg.snp_list) < max(1, conf.rlp_min_n_snps):
                continue
            if reg.snp_list[-1].pos - reg.snp_list[0].pos + 1 < conf.rlp_min_gap:
                continue
            if conf.debug > 2:
                debug("region '%s': do local phasing ..." % reg.name)
            dat = adata[:, (adata.var["chrom"] == reg.chrom) & \
                            (adata.var["pos"] >= reg.start) & \
                            (adata.var["pos"] < reg.end)].copy()
            reg, flip = reg_local_phasing(
                reg = reg,
                AD = dat.layers['AD'].copy(),
                DP = dat.layers['DP'].copy(),
                kws_localphase = None,
                verbose = True if conf.debug > 2 else False
            )
            if flip is None:
                n_rlp_failed += 1
                if conf.debug > 1:
                    debug("region '%s': local phasing failed." % reg.name)
            else:
                if conf.debug > 1:
                    debug("region '%s': #SNPs - total=%d; flipped=%d" % \
                          (reg.name, len(reg.snp_list), np.sum(flip)))
                n_slp_flipped += np.sum(flip)
            n_slp += len(reg.snp_list)
            n_rlp += 1
        del conf.snp_adata
        del conf.ref_cells
        gc.collect()
    info("#regions: total=%d; local_phasing=%d; local_phasing_failed=%d." % \
        (len(conf.reg_list), n_rlp, n_rlp_failed))
    info("#SNPs: local_phasing=%d; local_phasing_flipped=%d." % \
        (n_slp, n_slp_flipped))
    

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
    conf.snp_set.destroy()
    conf.snp_set = None

    
    # multiprocessing, push regions into process pool.
    thdata_list = []
    pool = multiprocessing.Pool(processes = m_thread)
    mp_result = []
    for i in range(m_thread):
        thdata = ThreadData(
            idx = i, conf = conf,
            reg_obj = reg_fn_list[i], is_reg_pickle = True,
            out_region_fn = conf.out_region_fn + "." + str(i),
            out_ad_fn = conf.out_ad_fn + "." + str(i),
            out_dp_fn = conf.out_dp_fn + "." + str(i),
            out_oth_fn = conf.out_oth_fn + "." + str(i),
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
        [td.out_ad_fn for td in thdata_list], ZF_F_GZIP, 
        conf.out_ad_fn, "w", ZF_F_PLAIN,
        nr_reg_list, len(conf.samples),
        sum([td.nr_ad for td in thdata_list]),
        remove = True
    ) < 0:
        raise ValueError("errcode -17")

    if merge_mtx(
        [td.out_dp_fn for td in thdata_list], ZF_F_GZIP, 
        conf.out_dp_fn, "w", ZF_F_PLAIN,
        nr_reg_list, len(conf.samples), 
        sum([td.nr_dp for td in thdata_list]),
        remove = True
    ) < 0:
        raise ValueError("errcode -19")

    if merge_mtx(
        [td.out_oth_fn for td in thdata_list], ZF_F_GZIP, 
        conf.out_oth_fn, "w", ZF_F_PLAIN,
        nr_reg_list, len(conf.samples),
        sum([td.nr_oth for td in thdata_list]),
        remove = True
    ) < 0:
        raise ValueError("errcode -21")



def afc_run(conf):
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
        ret = afc_core(conf)
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
        conf.out_dir, conf.out_prefix + "region.tsv")
    conf.out_sample_fn = os.path.join(
        conf.out_dir, conf.out_prefix + "samples.tsv")
    conf.out_ad_fn = os.path.join(conf.out_dir, conf.out_prefix + "AD.mtx")
    conf.out_dp_fn = os.path.join(conf.out_dir, conf.out_prefix + "DP.mtx")
    conf.out_oth_fn = os.path.join(conf.out_dir, conf.out_prefix + "OTH.mtx")

    
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

    
    if conf.snp_fn:
        if os.path.isfile(conf.snp_fn):
            if conf.snp_fn.endswith(".vcf") or conf.snp_fn.endswith(".vcf.gz")\
                    or conf.snp_fn.endswith(".vcf.bgz"):
                conf.snp_set = load_snp_from_vcf(conf.snp_fn, verbose = True)
            else:
                conf.snp_set = load_snp_from_tsv(conf.snp_fn, verbose = True)
            if not conf.snp_set or conf.snp_set.get_n() <= 0:
                error("failed to load snp file.")
                return(-1)
            else:
                info("%d SNPs loaded." % conf.snp_set.get_n())       
        else:
            error("snp file '%s' does not exist." % conf.snp_fn)
            return(-1)      
    else:
        error("SNP file needed!")
        return(-1)
    
    
    if conf.cellsnp_dir is not None:
        assert_e(conf.cellsnp_dir)
        snp_adata = csp_load_data(conf.cellsnp_dir)
        info("cellsnp SNP adata shape = %s." % str(snp_adata.shape))
        
        assert len(conf.samples) == snp_adata.shape[0]
        for cell in conf.samples:
            assert cell in snp_adata.obs['cell'].to_numpy()
        
        if snp_adata.shape[1] != conf.snp_set.get_n():
            warn("n_snp: snp_adata=%d; snp_set=%d!" % \
                (snp_adata.shape[1], conf.snp_set.get_n()))
            assert snp_adata.shape[1] >= conf.snp_set.get_n()
            
        idx_lst = []
        for i in range(snp_adata.shape[1]):
            chrom = snp_adata.var['chrom'].iloc[i]
            pos = snp_adata.var['pos'].iloc[i]
            hits = conf.snp_set.fetch(
                chrom = chrom,
                start = pos,
                end = pos + 1
            )
            if len(hits) > 0:     # in case SNPs were filtered in phasing.
                idx_lst.append(i)
            else:
                if conf.debug > 0:
                    warn("SNP '%s:%d' was filtered before!" % (chrom, pos))
                
        if len(idx_lst) < snp_adata.shape[1]:
            try:
                snp_adata = snp_adata[:, snp_adata.var.index.iloc[idx_lst]].copy()
            except:
                snp_adata = snp_adata[:, snp_adata.var.index.take(idx_lst)].copy()
            info("SNP adata shape after subset: %s." % str(snp_adata.shape))
        conf.snp_adata = snp_adata.copy()


    if conf.ref_cell_fn is not None:
        assert_e(conf.ref_cell_fn)
        conf.ref_cells = np.genfromtxt(
            conf.ref_cell_fn, dtype = "str", delimiter = "\t")
        assert len(conf.ref_cells) <= len(conf.samples)
        for cell in conf.ref_cells:
            assert cell in conf.samples
        
        
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
