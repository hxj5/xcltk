# zfile.py - wrapper of file object, supporting GZIP/BGZIP. 

# Note that there's a bug when reading BGZIP file with pysam (issue #438):
#   it will stop reading when encounting a blank line, this could be caused 
#   by the issue that `bgzf_getline()` function always drops the tailing '\n'.
#   This bug exists in pysam v0.18.0 and before(?).
# So we need to use gzip to read GZIP and BGZIP files
# it's probably ok to use pysam to write BGZIP files.

import gzip
import pysam


class ZFile:
    """Simple File object wrapper that supports plain/gzip/bgzip.

    Attributes
    ----------
    file_name : str
        File name.
    mode : str
        File mode.
    file_type : int
        File type / format, one of `ZF_F_XXX`.
    is_bytes : bool
        Whether the data is of type `bytes`.
    encoding : str
        Encoding for `bytes` data; set to "utf8" if None.
    """
    def __init__(self, file_name, mode, file_type, 
                 is_bytes = False, encoding = None):
        self.file_name = file_name
        self.file_type = file_type
        self.mode = mode
        self.is_bytes = is_bytes
        self.encoding = encoding if encoding else "utf8"

        self.fp = None
        self.buf = self.__reset_buf()

        if file_type == ZF_F_AUTO:
            fn = file_name.lower()
            if fn.endswith(".gz") or fn.endswith(".gzip"):
                file_type = ZF_F_GZIP
            else:
                file_type = ZF_F_PLAIN

        if file_type == ZF_F_PLAIN:
            self.fp = open(self.file_name, self.mode)
        elif file_type == ZF_F_GZIP:
            self.fp = gzip.open(self.file_name, self.mode)
        elif file_type == ZF_F_BGZIP:
            if "r" in self.mode:
                self.fp = gzip.open(self.file_name, self.mode)
            else:
                self.fp = pysam.BGZFile(self.file_name, self.mode)
        else:
            raise ValueError("invalid file type")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __iter__(self):
        return self.fp

    def __next__(self):
        line = self.readline()
        if not line:
            raise StopIteration()
        return(line)

    def __reset_buf(self):
        return(bytes("", self.encoding) if self.is_bytes else "")

    def close(self):
        if self.fp:
            if self.buf:
                self.fp.write(self.buf)
                self.buf = None
            self.fp.close()
            self.fp = None

    def read(self, size = None):
        if not self.fp:
            raise OSError()
        if size is None:
            return(self.fp.read())
        else:
            return(self.fp.read(size))

    def readline(self, size = None):
        if not self.fp:
            raise OSError()
        if size is None:
            return(self.fp.readline())
        else:
            return(self.fp.readline(size))

    def readlines(self, size = None):
        if not self.fp:
            raise OSError()
        if size is None:
            return(self.fp.readlines())
        else:
            return(self.fp.readlines(size))

    def write(self, data):
        if not self.fp:
            raise OSError()
        self.buf += data
        if len(self.buf) >= ZF_BUFSIZE:
            ret = self.fp.write(self.buf)
            self.buf = self.__reset_buf()
            return(ret)
        return(len(data))

def zopen(file_name, mode, file_type = None, is_bytes = False, encoding = None):
    if not file_name:
        raise OSError()
    if file_type is None:
        file_type = ZF_F_AUTO
    return ZFile(file_name, mode, file_type, is_bytes, encoding)

ZF_F_PLAIN = 0
ZF_F_GZIP = 1
ZF_F_BGZIP = 2
ZF_F_AUTO = 3

ZF_BUFSIZE = 1048576   # 1M

# TODO: update the debugging codes below
if __name__ == "__main__":
    pass
