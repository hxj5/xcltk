# thread.py - wrapper of data for each thread.

import sys

class ThreadData:
    """Thread Data"""
    def __init__(self, 
        idx, conf, 
        reg_obj, is_reg_pickle, 
        out_region_fn,
        out_mtx_fn,
        out_fn = None
    ):
        self.idx = idx
        self.conf = conf

        self.reg_obj = reg_obj
        self.is_reg_pickle = is_reg_pickle

        self.out_region_fn = out_region_fn
        self.out_mtx_fn = out_mtx_fn

        self.out_fn = out_fn

        self.nr_reg = 0
        self.nr_mtx = 0
        
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
        s += "%sout_mtx_fn = %s\n" % (prefix, self.out_mtx_fn)

        s += "%sout_fn = %s\n" % (prefix, self.out_fn)
        
        s += "%snum_record_region = %d\n" % (prefix, self.nr_reg)
        s += "%snum_record_mtx = %d\n" % (prefix, self.nr_mtx)

        s += "%sreturn code = %d\n" % (prefix, self.ret)
        s += "%s\n" % prefix

        fp.write(s)