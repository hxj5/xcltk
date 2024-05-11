# thread.py - wrapper of data for each thread.

import sys

class ThreadData:
    """Thread Data"""
    def __init__(self, 
        idx, conf, 
        reg_obj, is_reg_pickle, 
        out_region_fn,
        out_ad_fn, out_dp_fn, out_oth_fn,
        out_fn = None
    ):
        self.idx = idx
        self.conf = conf

        self.reg_obj = reg_obj
        self.is_reg_pickle = is_reg_pickle

        self.out_region_fn = out_region_fn
        self.out_ad_fn = out_ad_fn
        self.out_dp_fn = out_dp_fn
        self.out_oth_fn = out_oth_fn

        self.out_fn = out_fn

        self.nr_reg = 0
        self.nr_ad = 0
        self.nr_dp = 0
        self.nr_oth = 0
        
        self.ret = -1

    def destroy(self):
        pass

    def show(self, fp = None, prefix = ""):
        if fp is None:
            fp = sys.stderr
        
        s = "%s\n" % prefix
        s += "%sindex = %d\n" % (prefix, self.idx)

        s += "%sis_reg_pickle = %s\n" % (prefix, self.is_reg_pickle)
        if self.is_reg_pickle:
            s += "%sreg_obj filename = %s\n" % (prefix, self.reg_obj)
        else:
            s += "%slen(reg_obj) = %d\n" % (prefix, len(self.reg_obj))

        s += "%sout_region_fn = %s\n" % (prefix, self.out_region_fn)
        s += "%sout_ad_fn = %s\n" % (prefix, self.out_ad_fn)
        s += "%sout_dp_fn = %s\n" % (prefix, self.out_dp_fn)
        s += "%sout_oth_fn = %s\n" % (prefix, self.out_oth_fn)

        s += "%sout_fn = %s\n" % (prefix, self.out_fn)
        
        s += "%snum_record_region = %d\n" % (prefix, self.nr_reg)
        s += "%snum_record_ad = %d\n" % (prefix, self.nr_ad)
        s += "%snum_record_dp = %d\n" % (prefix, self.nr_dp)
        s += "%snum_record_oth = %d\n" % (prefix, self.nr_oth)

        s += "%sreturn code = %d\n" % (prefix, self.ret)
        s += "%s\n" % prefix

        fp.write(s)
