# Base Utils
# Author: Xianjie Huang

import os
import sys
import time

def assert_e(path, name, _type = None):
    """
    @abstract     Assert path exists and exit the program if not.
                  Inspired by test -e in shell
    @param path   Path to the file/dir [str]
    @param name   Path name [str]
    @param _type  Path type: file or dir. if None, auto detect [str]
    @return       Void
    @example      assert_e("test.sh", "file", "test script")
    """
    if not path:
        raise ValueError("path is empty for '%s'" % name)
    exist = False
    if _type: 
        if _type.lower() == "file":
            if os.path.isfile(path):
                exist = True
        elif _type.lower() == "dir":
            if os.path.isdir(path):
                exist = True
        else:
            raise ValueError("path type should be file or dir for '%s'" % name)
    else:
        if os.path.isfile(path) or os.path.isdir(path):
            exist = True
    if not exist:
        raise IOError("%s '%s' does not exist." % (name, path))

def assert_n(var, name):
    """
    @abstract    Assert variable is valid, i.e., True for`if <variable>`.
                 Exit the program if not. Inspired by test -n in shell
    @param var   The variable
    @param name  Variable name [str]
    @return      Void
    @example     assert_n("hello", "test variable")
    """
    if not var:
        raise ValueError("var is empty for '%s'" % name)

def get_now_str(fmt = "%Y-%m-%d %H:%M:%S"):
    """
    @abstract   Return string of now
    @param fmt  Time string format [str]
    @return     String of now [str] 
    """
    return time.strftime(fmt, time.localtime())

def log(msg, fp = sys.stdout):
    """
    @abstract   Format log message and print
    @param msg  Log message to be printed [str]
    @param fp   File pointer [FILE*]
    @return     Void
    """
    fp.write("[%s] %s\n" % (get_now_str(), msg))