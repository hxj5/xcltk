# xcltk.py - xcltk cmdline interface.

import sys

from .baf.count import afc_main as baf_allelefc
from .baf.fixref import fixref_main as baf_fixref
from .baf.pipeline import pipeline_main as baf_baf
from .baf.rpc import rpc_main as baf_rpc
from .config import APP, VERSION
from .rdr.count import fc_main as rdr_basefc
from .tools.convert import convert_main as tools_convert


def __usage(fp = sys.stdout):
    msg =  "\n"
    msg += "Program: %s (Toolkit for XClone Preprocessing)\n" % APP
    msg += "Version: %s\n" % VERSION
    msg += "\n"
    msg += "Usage:   %s <command> [options]\n" % APP
    msg += "\n"                                          \
           "Commands:\n"                                  \
           "  -- BAF calculation\n"                                         \
           "     allelefc         Allele-specific feature counting.\n"      \
           "     baf              Preprocessing pipeline for XClone BAF.\n"   \
           "     fixref           Fix REF allele mismatches based on reference FASTA.\n"  \
           "     rpc              Reference phasing correction.\n"            \
           "\n"                                                              \
           "  -- RDR calculation\n"                                          \
           "     basefc           Basic feature counting.\n"                 \
           "\n"                                                              \
           "  -- Tools\n"                                        \
           "     convert          Convert between different formats of genomic features.\n"   \
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
        sys.exit(0)

    command = sys.argv[1]
    if command == "allelefc": baf_allelefc(sys.argv)
    elif command == "baf": baf_baf(sys.argv)
    elif command == "fixref": baf_fixref(sys.argv)
    elif command == "rpc": baf_rpc(sys.argv)
    elif command == "basefc": rdr_basefc(sys.argv)
    elif command == "convert": tools_convert(sys.argv)
    elif command in ("-h", "--help"): __usage(); sys.exit(0)
    elif command in ("-V", "--version"): sys.stderr.write("%s\n" % VERSION); sys.exit(0)
    else: sys.stderr.write("Error: wrong command '%s'\n" % command); sys.exit(1)