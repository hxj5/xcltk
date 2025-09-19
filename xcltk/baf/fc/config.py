# config.py - configuration


import sys
from ...config import APP


class Config:
    def __init__(self):
        self.defaults = DefaultConfig()
        self.argv = None

        self.sam_fn = None
        self.sam_list_fn = None
        self.barcode_fn = None
        self.sample_id_str = None
        self.sample_id_fn = None
        self.region_fn = None
        self.snp_fn = None
        self.out_dir = None
        self.debug = self.defaults.DEBUG

        self.cellsnp_dir = None
        self.ref_cell_fn = None
        self.cell_tag = self.defaults.CELL_TAG
        self.umi_tag = self.defaults.UMI_TAG
        self.nproc = self.defaults.NPROC
        self.min_count = self.defaults.MIN_COUNT
        self.min_maf = self.defaults.MIN_MAF
        self.output_all_reg = self.defaults.OUTPUT_ALL_REG
        self.no_dup_hap = self.defaults.NO_DUP_HAP

        self.min_mapq = self.defaults.MIN_MAPQ
        self.min_len = self.defaults.MIN_LEN
        self.incl_flag = self.defaults.INCL_FLAG
        self.excl_flag = -1
        self.no_orphan = self.defaults.NO_ORPHAN


        # derived parameters.
        
        self.barcodes = None     # list of barcode strings.
        self.sample_ids = None
        self.reg_list = None     # list of gene/block objects.
        self.snp_set = None      # set of SNPs.

        self.sam_fn_list = None
        self.samples = None
        
        self.snp_adata = None
        self.ref_cells = None


        # internal parameters
        
        self.out_prefix = APP + "."
        self.out_region_fn = None
        self.out_sample_fn = None
        self.out_ad_fn = None
        self.out_dp_fn = None
        self.out_oth_fn = None
        
        # rlp: region wise local phasing
        # - requirements below must be meeted before doing local phasing in 
        #   one region.
        self.rlp_min_len = 50000       # minimum length of the region.
        self.rlp_min_n_snps = 2        # minimum number of expressed SNPs.
        self.rlp_min_gap = 50000       # minimum gap between the first SNP and last SNP.

        
    def show(self, fp = None, prefix = ""):
        if fp is None:
            fp = sys.stderr

        s =  "%s\n" % prefix
        s += "%ssam_file = %s\n" % (prefix, self.sam_fn)
        s += "%ssam_list_file = %s\n" % (prefix, self.sam_list_fn)
        s += "%sbarcode_file = %s\n" % (prefix, self.barcode_fn)
        s += "%ssample_id_str = %s\n" % (prefix, self.sample_id_str)
        s += "%ssample_id_file = %s\n" % (prefix, self.sample_id_fn)
        s += "%sregion_file = %s\n" % (prefix, self.region_fn)
        s += "%ssnp_file = %s\n" % (prefix, self.snp_fn)
        s += "%sout_dir = %s\n" % (prefix, self.out_dir)
        s += "%sdebug_level = %d\n" % (prefix, self.debug)
        s += "%s\n" % prefix

        s += "%scellsnp_dir = %s\n" % (prefix, self.cellsnp_dir)
        s += "%sref_cell_fn = %s\n" % (prefix, self.ref_cell_fn)
        s += "%scell_tag = %s\n" % (prefix, self.cell_tag)
        s += "%sumi_tag = %s\n" % (prefix, self.umi_tag)
        s += "%snumber_of_processes = %d\n" % (prefix, self.nproc)
        s += "%smin_count = %d\n" % (prefix, self.min_count)
        s += "%smin_maf = %f\n" % (prefix, self.min_maf)
        s += "%soutput_all_reg = %s\n" % (prefix, self.output_all_reg)
        s += "%sno_dup_hap = %s\n" % (prefix, self.no_dup_hap)
        s += "%s\n" % prefix

        s += "%smin_mapq = %d\n" % (prefix, self.min_mapq)
        s += "%smin_len = %d\n" % (prefix, self.min_len)
        s += "%sinclude_flag = %d\n" % (prefix, self.incl_flag)
        s += "%sexclude_flag = %d\n" % (prefix, self.excl_flag)
        s += "%sno_orphan = %s\n" % (prefix, self.no_orphan)
        s += "%s\n" % prefix

        s += "%s#BAMs = %d\n" % (prefix, len(self.sam_fn_list) if \
                self.sam_fn_list is not None else -1)
        s += "%s#barcodes = %d\n" % (prefix, len(self.barcodes) if \
                self.barcodes is not None else -1)
        s += "%s#sample IDs = %d\n" % (prefix, len(self.sample_ids) \
                if self.sample_ids is not None else -1)
        s += "%s#regions = %d\n" % (prefix, len(self.reg_list) if \
                self.reg_list is not None else -1)
        s += "%s#snps = %d\n" % (prefix, self.snp_set.get_n() if \
                self.snp_set is not None else -1)
        s += "%sshape of snp adata = %s\n" % (prefix, str(self.snp_adata.shape)\
                if self.snp_adata is not None else 'None')
        s += "%s#reference cells = %d\n" % (prefix, len(self.ref_cells) \
                if self.ref_cells is not None else -1)
        s += "%s\n" % prefix

        s += "%soutput_region_file = %s\n" % (prefix, self.out_region_fn)
        s += "%soutput_sample_file = %s\n" % (prefix, self.out_sample_fn)
        s += "%soutput_ad_file = %s\n" % (prefix, self.out_ad_fn)
        s += "%soutput_dp_file = %s\n" % (prefix, self.out_dp_fn)
        s += "%soutput_oth_file = %s\n" % (prefix, self.out_oth_fn)
        s += "%s\n" % prefix
        
        s += "%srlp_min_len = %d\n" % (prefix, self.rlp_min_len)
        s += "%srlp_min_n_snps = %d\n" % (prefix, self.rlp_min_n_snps)
        s += "%srlp_min_gap = %d\n" % (prefix, self.rlp_min_gap)
        s += "%s\n" % prefix

        fp.write(s)

        
    def use_barcodes(self):
        return self.cell_tag is not None
    
    def use_local_phasing(self):
        return self.cellsnp_dir is not None

    def use_umi(self):
        return self.umi_tag is not None



class DefaultConfig:
    def __init__(self):
        self.DEBUG = 0
        self.CELL_TAG = "CB"
        self.UMI_TAG = "UB"
        self.UMI_TAG_BC = "UB"    # the default umi tag for 10x data.
        self.NPROC = 1
        self.MIN_COUNT = 1 
        self.MIN_MAF = 0
        self.OUTPUT_ALL_REG = False
        self.NO_DUP_HAP = True

        self.MIN_MAPQ = 20
        self.MIN_LEN = 30
        self.INCL_FLAG = 0
        self.EXCL_FLAG_UMI = 772
        self.EXCL_FLAG_XUMI = 1796
        self.NO_ORPHAN = True



if __name__ == "__main__":
    conf = Config()
    conf.show()
