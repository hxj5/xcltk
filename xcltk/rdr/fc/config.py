# config.py - configuration.

import sys

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
        self.out_dir = None
        self.debug = self.defaults.DEBUG

        self.cell_tag = self.defaults.CELL_TAG
        self.umi_tag = self.defaults.UMI_TAG
        self.nproc = self.defaults.NPROC
        self.output_all_reg = self.defaults.OUTPUT_ALL_REG

        self.min_mapq = self.defaults.MIN_MAPQ
        self.min_len = self.defaults.MIN_LEN
        self.min_include = self.defaults.MIN_INCLUDE
        self.incl_flag = self.defaults.INCL_FLAG
        self.excl_flag = -1
        self.no_orphan = self.defaults.NO_ORPHAN

        self.barcodes = None     # list of barcode strings.
        self.sample_ids = None
        self.reg_list = None     # list of gene/block objects.

        self.sam_fn_list = None
        self.samples = None

        self.out_prefix = ""
        self.out_region_fn = None
        self.out_sample_fn = None
        self.out_mtx_fn = None

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
        s += "%sout_dir = %s\n" % (prefix, self.out_dir)
        s += "%sdebug_level = %d\n" % (prefix, self.debug)
        s += "%s\n" % prefix

        s += "%scell_tag = %s\n" % (prefix, self.cell_tag)
        s += "%sumi_tag = %s\n" % (prefix, self.umi_tag)
        s += "%snumber_of_processes = %d\n" % (prefix, self.nproc)
        s += "%soutput_all_reg = %s\n" % (prefix, self.output_all_reg)
        s += "%s\n" % prefix

        s += "%smin_mapq = %d\n" % (prefix, self.min_mapq)
        s += "%smin_len = %d\n" % (prefix, self.min_len)
        s += "%smin_include = %f\n" % (prefix, self.min_include)
        s += "%sinclude_flag = %d\n" % (prefix, self.incl_flag)
        s += "%sexclude_flag = %d\n" % (prefix, self.excl_flag)
        s += "%sno_orphan = %s\n" % (prefix, self.no_orphan)
        s += "%s\n" % prefix

        s += "%snumber_of_BAMs = %d\n" % (prefix, len(self.sam_fn_list) if \
                self.sam_fn_list is not None else -1)
        s += "%snumber_of_barcodes = %d\n" % (prefix, len(self.barcodes) if \
                self.barcodes is not None else -1)
        s += "%snumber_of_sample_IDs = %d\n" % (prefix, len(self.sample_ids) \
                if self.sample_ids is not None else -1)
        s += "%snumber_of_regions = %d\n" % (prefix, len(self.reg_list) if \
                self.reg_list is not None else -1)
        s += "%s\n" % prefix

        s += "%soutput_region_file = %s\n" % (prefix, self.out_region_fn)
        s += "%soutput_sample_file = %s\n" % (prefix, self.out_sample_fn)
        s += "%soutput_mtx_file = %s\n" % (prefix, self.out_mtx_fn)
        s += "%s\n" % prefix

        fp.write(s)

    def use_barcodes(self):
        return self.cell_tag is not None

    def use_umi(self):
        return self.umi_tag is not None


class DefaultConfig:
    def __init__(self):
        self.DEBUG = 0
        self.CELL_TAG = "CB"
        self.UMI_TAG = "UB"
        self.UMI_TAG_BC = "UB"    # the default umi tag for 10x data.
        self.NPROC = 1
        self.OUTPUT_ALL_REG = True

        self.MIN_MAPQ = 20
        self.MIN_LEN = 30
        self.MIN_INCLUDE = 0.9
        self.INCL_FLAG = 0
        self.EXCL_FLAG_UMI = 772
        self.EXCL_FLAG_XUMI = 1796
        self.NO_ORPHAN = True


if __name__ == "__main__":
    conf = Config()
    conf.show()
