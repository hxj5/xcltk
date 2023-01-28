# cmdline options
# Author: Xianjie Huang 

#TODO: 
# - basefc: support multiple bam files
# - add install list (conda + pip)

import sys
from .config import PROGRAM, VERSION
from .baf.fixref import fixref as baf_fix_ref
from .baf.pileup import pileup as baf_pileup
from .baf.rpc import main_rpc as baf_rpc
from .rdr.basefc import base_fc as rdr_base_fc
from .region.convert import convert as reg_convert

def __usage(fp = sys.stderr):
    msg =  "\n"
    msg += "Program: %s (Toolkit for XClone)\n" % PROGRAM
    msg += "Version: %s\n" % VERSION
    msg += "\n"
    msg += "Usage:   %s <command> [options]\n" % PROGRAM
    msg += "\n"                                          \
           "Commands:\n"                                  \
           "  -- BAF calculation\n"                                         \
           "     fixref           Fix REF, ALT and GT.\n"                    \
           "     rpc              Reference phasing correction.\n"            \
           "     pileup           Pileup, support unique counting.\n"         \
           "\n"                                                              \
           "  -- RDR calculation\n"                                          \
           "     basefc           Basic feature counting.\n"                 \
           "\n"                                                              \
           "  -- Region operations\n"                                        \
           "     convert          Convert different region file formats.\n"   \
           "\n"                                                              \
           "  -- Others\n"                                                   \
           "     -h, --help       Print this message and exit.\n"             \
           "     -V, --version    Print version and exit.\n"                 \
           "\n"
    fp.write(msg)

def main():
    argc = len(sys.argv)
    if argc < 2:
        __usage()
        sys.exit(1)

    command = sys.argv[1]
    if command == "fixref": baf_fix_ref(sys.argv)
    elif command == "rpc": baf_rpc(sys.argv)
    elif command == "pileup": baf_pileup(sys.argv)
    elif command == "basefc": rdr_base_fc(sys.argv)
    elif command == "convert": reg_convert(sys.argv)
    elif command in ("-h", "--help"): __usage(); sys.exit(3)
    elif command in ("-V", "--version"): sys.stderr.write("%s\n" % VERSION); sys.exit(3)
    else: sys.stderr.write("Error: wrong command '%s'\n" % command); sys.exit(5)

if __name__ == "__main__":
    main()

