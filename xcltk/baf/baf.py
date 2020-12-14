# cmdline options
# Author: Xianjie Huang 

import sys
from .phase_snp import phase_snp
from .config import APP

def __usage(fp = sys.stderr):
    msg =  "\n"
    msg += "Usage: %s <command> [options]\n" % APP
    msg += "\n"                                          \
           "Commands:\n"                                  \
           "  phase_snp      Aggregate SNPs into haplotype blocks.\n"     \
           "  -h, --help     Print this message.\n"                      \
           "\n"
    fp.write(msg)

def main():
    argc = len(sys.argv)
    if argc < 2:
        __usage()
        sys.exit(1)

    command = sys.argv[1]
    if command == "phase_snp":  phase_snp(argc, sys.argv)
    elif command in ("-h", "--help"): __usage(); sys.exit(3)
    else: sys.stderr.write("Error: wrong command '%s'\n" % command); sys.exit(5)

if __name__ == "__main__":
    main()

