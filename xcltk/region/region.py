# cmdline options
# Author: Xianjie Huang 

import sys
from .convert import convert
from .config import APP

def __usage(fp = sys.stderr):
    msg =  "\n"
    msg += "Usage: %s <command> [options]\n" % APP 
    msg += "\n"                                          \
           "Commands:\n"                                  \
           "  convert        Convert different region file formats.\n"  \
           "  -h, --help     Print this message.\n"                    \
           "\n"
    fp.write(msg)

def main():
    argc = len(sys.argv)
    if argc < 2:
        __usage()
        sys.exit(1)

    command = sys.argv[1]
    if command == "convert": convert()
    elif command in ("-h", "--help"): __usage(); sys.exit(3)
    else: sys.stderr.write("Error: wrong command '%s'\n" % command); sys.exit(5)

if __name__ == "__main__":
    main()

