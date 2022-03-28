# Global Configure

import sys
from .app import APP

class Config:
    def __init__(self):
        self.sam_fn = None
        self.sam_list_fn = None
        self.out_dir = None
        self.region_fn = None
        self.barcode_fn = None
        self.sid = None
        self.sid_fn = None
        self.debug = CFG_DEBUG

        self.cell_tag = CFG_CELL_TAG
        self.umi_tag = CFG_UMI_TAG
        self.nproc = CFG_NPROC
        self.min_count = CFG_MIN_COUNT
        self.pair_end = CFG_PAIR_END

        self.min_overlap = CFG_MIN_OVERLAP
        self.max_gap = CFG_MAX_GAP
        self.min_uread = CFG_MIN_UREAD
        self.min_ufrac = CFG_MIN_UFRAC
        self.min_jumi = CFG_MIN_JUMI
        self.min_jcell = CFG_MIN_JCELL
        self.incl_denovo = CFG_INCL_DENOVO
        self.denovo_jc_len = CFG_DENOVO_JC_LEN

        self.min_mapq = CFG_MIN_MAPQ
        self.min_len = CFG_MIN_LEN
        self.incl_flag = CFG_INCL_FLAG
        self.excl_flag = CFG_EXCL_FLAG_UMI
        self.no_orphan = CFG_NO_ORPHAN

        self.sam_list = None     # list of pysam::AlignmentFile objects.
        self.sam_fn_list = None  # list of sam filenames.
        self.barcodes = None     # list of barcode strings.
        self.sid_list = None     # list of sample ID strings.
        self.reg_list = None     # list of gene objects.

        self.out_prefix = APP + "."
        self.out_junction_fn = None
        self.out_region_fn = None
        self.out_sample_fn = None
        self.out_spliced_fn = None
        self.out_unspliced_fn = None
        self.out_ambiguous_fn = None

    def destroy(self):
        if self.sam_list:
            for sam in self.sam_list:
                sam.close()
            self.sam_list.clear()
        if self.sam_fn_list:
            self.sam_fn_list.clear()
        if self.barcodes:
            self.barcodes.clear()
        if self.sid_list:
            self.sid_list.clear()
        if self.reg_list:
            for reg in self.reg_list:
                reg.destroy()
            self.reg_list.clear()

    def show(self, fp = None, prefix = ""):
        if fp is None:
            fp = sys.stderr

        s =  "%s\n" % prefix
        s += "%ssam_file = %s\n" % (prefix, self.sam_fn)
        s += "%ssam_list_file = %s\n" % (prefix, self.sam_list_fn)
        s += "%sout_dir = %s\n" % (prefix, self.out_dir)
        s += "%sregion_file = %s\n" % (prefix, self.region_fn)
        s += "%sbarcode_file = %s\n" % (prefix, self.barcode_fn)
        s += "%ssample_id (str) = %s\n" % (prefix, self.sid)
        s += "%ssample_id_list_file = %s\n" % (prefix, self.sid_fn)
        s += "%sdebug = %d\n" % (prefix, self.debug)
        s += "%s\n" % prefix

        s += "%scell_tag = %s\n" % (prefix, self.cell_tag)
        s += "%sumi_tag = %s\n" % (prefix, self.umi_tag)
        s += "%snumber_of_processes = %d\n" % (prefix, self.nproc)
        s += "%smin_count = %d\n" % (prefix, self.min_count)
        s += "%spair_end = %s\n" % (prefix, self.pair_end)

        s += "%smin_overlap = %d\n" % (prefix, self.min_overlap)
        s += "%smax_gap = %d\n" % (prefix, self.max_gap)
        s += "%smin_uread = %d\n" % (prefix, self.min_uread)
        s += "%smin_ufrac = %f\n" % (prefix, self.min_ufrac)
        s += "%smin_jumi = %d\n" % (prefix, self.min_jumi)
        s += "%smin_jcell = %d\n" % (prefix, self.min_jcell)
        s += "%sincl_denovo = %s\n" % (prefix, self.incl_denovo)
        s += "%sdenovo_jc_len = %f\n" % (prefix, self.denovo_jc_len)
        s += "%s\n" % prefix

        s += "%smin_mapq = %d\n" % (prefix, self.min_mapq)
        s += "%smin_len = %d\n" % (prefix, self.min_len)
        s += "%sinclude_flag = %d\n" % (prefix, self.incl_flag)
        s += "%sexclude_flag = %d\n" % (prefix, self.excl_flag)
        s += "%sno_orphan = %s\n" % (prefix, self.no_orphan)
        s += "%s\n" % prefix

        s += "%snumber_of_sams = %d\n" % (prefix, len(self.sam_fn_list) if \
                self.sam_fn_list is not None else -1)
        s += "%snumber_of_barcodes = %d\n" % (prefix, len(self.barcodes) if \
                self.barcodes is not None else -1)
        s += "%snumber_of_sample_IDs = %d\n" % (prefix, len(self.sid_list) if \
                self.sid_list is not None else -1)
        s += "%snumber_of_regions = %d\n" % (prefix, len(self.reg_list) if \
                self.reg_list is not None else -1)
        s += "%s\n" % prefix

        s += "%soutput_junction_file = %s\n" % (prefix, self.out_junction_fn)
        s += "%soutput_region_file = %s\n" % (prefix, self.out_region_fn)
        s += "%soutput_sample_file = %s\n" % (prefix, self.out_sample_fn)
        s += "%soutput_spliced_file = %s\n" % (prefix, self.out_spliced_fn)
        s += "%soutput_unspliced_file = %s\n" % (prefix, self.out_unspliced_fn)
        s += "%soutput_ambiguous_file = %s\n" % (prefix, self.out_ambiguous_fn)
        s += "%s\n" % prefix

        fp.write(s)

    def use_barcode(self):
        return self.cell_tag and self.barcodes

    def use_umi(self):
        return self.umi_tag is not None

CFG_DEBUG = 0
CFG_CELL_TAG = "CB"
#<DEV/>#CFG_UMI_TAG = "Auto"
CFG_UMI_TAG = "UB"
CFG_UMI_TAG_BC = "UB"    # the default umi tag for 10x data.
CFG_NPROC = 1
CFG_MIN_COUNT = 1        # not used yet
CFG_PAIR_END = False

CFG_MIN_OVERLAP = 5
CFG_MAX_GAP = 1
CFG_MIN_UREAD = 2
CFG_MIN_UFRAC = 2.0 / 3
CFG_MIN_JUMI = 2
CFG_MIN_JCELL = 1
CFG_INCL_DENOVO = False
CFG_DENOVO_JC_LEN = 0

CFG_MIN_MAPQ = 20
CFG_MIN_LEN = 30
CFG_INCL_FLAG = 0
CFG_EXCL_FLAG_UMI = 772
CFG_EXCL_FLAG_XUMI = 1796
CFG_NO_ORPHAN = True

if __name__ == "__main__":
    conf = Config()
    conf.show()
