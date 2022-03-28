# Thread
# Author: Xianjie Huang

import sys

class ThreadData:
    """Thread Data"""
    def __init__(self, 
        idx, conf, 
        reg_obj, is_reg_pickle, 
        out_junction_fn, out_region_fn,
        out_spliced_fn, out_unspliced_fn, out_ambiguous_fn,
        out_fn = None
    ):
        self.idx = idx
        self.conf = conf

        self.reg_obj = reg_obj
        self.is_reg_pickle = is_reg_pickle

        self.out_junction_fn = out_junction_fn
        self.out_region_fn = out_region_fn
        self.out_spliced_fn = out_spliced_fn
        self.out_unspliced_fn = out_unspliced_fn
        self.out_ambiguous_fn = out_ambiguous_fn

        self.out_fn = out_fn

        #self.nr_jc = 0
        self.nr_reg = 0
        self.nr_sp = 0
        self.nr_us = 0
        self.nr_am = 0
        
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

        s += "%sout_junction_fn = %s\n" % (prefix, self.out_junction_fn)
        s += "%sout_region_fn = %s\n" % (prefix, self.out_region_fn)
        s += "%sout_spliced_fn = %s\n" % (prefix, self.out_spliced_fn)
        s += "%sout_unspliced_fn = %s\n" % (prefix, self.out_unspliced_fn)
        s += "%sout_ambiguous_fn = %s\n" % (prefix, self.out_ambiguous_fn)

        s += "%sout_fn = %s\n" % (prefix, self.out_fn)
        
        #s += "%snum_record_junction_fn = %d\n" % (prefix, self.nr_jc)
        s += "%snum_record_region_fn = %d\n" % (prefix, self.nr_reg)
        s += "%snum_record_spliced_fn = %d\n" % (prefix, self.nr_sp)
        s += "%snum_record_unspliced_fn = %d\n" % (prefix, self.nr_us)
        s += "%snum_record_ambiguous_fn = %d\n" % (prefix, self.nr_am)

        s += "%sreturn code = %d\n" % (prefix, self.ret)
        s += "%s\n" % prefix

        fp.write(s)
