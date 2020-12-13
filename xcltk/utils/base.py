# Base Utils
# Author: Yuanhua Huang
# Modified by Xianjie Huang

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

def id_mapping(IDs1, IDs2, uniq_ref_only=True, IDs2_sorted=False):
    """
    Mapping IDs2 to IDs1. IDs1 (ref id) can have repeat values, but IDs2 need 
    to only contain unique ids.
    Therefore, IDs2[rv_idx] will be the same as IDs1.
    
    Parameters
    ----------
    IDs1 : array_like or list
        ids for reference.
    IDs2 : array_like or list
        ids waiting to map.
        
    Returns
    -------
    RV_idx : array_like, the same length of IDs1
        The index for IDs2 mapped to IDs1. If an id in IDs1 does not exist 
        in IDs2, then return a None for that id.
    """
    idx1 = sorted(range(len(IDs1)), key=IDs1.__getitem__)
    if IDs2_sorted:
        idx2 = range(len(IDs2))
    else:
        idx2 = sorted(range(len(IDs2)), key=IDs2.__getitem__)
    RV_idx1, RV_idx2 = [], []
    
    i, j = 0, 0
    while i < len(idx1):
        if j == len(idx2) or IDs1[idx1[i]] < IDs2[idx2[j]]:
            RV_idx1.append(idx1[i])
            RV_idx2.append(None)
            i += 1
        elif IDs1[idx1[i]] == IDs2[idx2[j]]:
            RV_idx1.append(idx1[i])
            RV_idx2.append(idx2[j])
            i += 1
            if uniq_ref_only: j += 1
        elif IDs1[idx1[i]] > IDs2[idx2[j]]:
            j += 1
            
    origin_idx = sorted(range(len(RV_idx1)), key=RV_idx1.__getitem__)
    RV_idx = [RV_idx2[i] for i in origin_idx]
    return RV_idx

def unique_list(X):
    """unique a list with index and count
    Example
    -------
    >>> unique_list([1,2,4,5,3,2,4])
    >>> ([1, 2, 3, 4, 5], [0, 1, 4, 2, 3], [1, 2, 1, 2, 1])
    """
    idx = sorted(range(len(X)), key=X.__getitem__)
    X_uniq = []
    X_count = []
    idx_uniq = []
    for i in idx:
        if len(X_uniq) == 0 or X[i] != X_uniq[-1]:
            X_uniq.append(X[i])
            X_count.append(1)
            idx_uniq.append(i)
        else:
            X_count[-1] += 1
    return X_uniq, idx_uniq, X_count

