# Basic Utils
# Author: Xianjie Huang <hxj5@hku.hk>

import os
import sys
import time

def assert_path_exists(path, _type, name):
    """
    @abstract     Assert path exists and exit the program if not
    @param path   Path to the file/dir [str]
    @param _type  Path type: file or dir [str]
    @param name   Path name [str]
    @return       Void
    @example      check_path_exists("test.sh", "file", "test script")
    """
    if not path:
        raise ValueError("path is empty")
    exist = False
    if _type and _type.lower() == "file":
        if os.path.isfile(path):
            exist = True
    elif _type and _type.lower() == "dir":
        if os.path.isdir(path):
            exist = True
    else:
        raise ValueError("path type should be file or dir.")
    if not exist:
        raise IOError("%s '%s' does not exist." % (name, path))

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

