# convert.py - convert between different formats of genomic features.

#TODO:
#- use getopt to re-implement cmdline.

import os
import sys

from optparse import OptionParser

from ..config import APP, VERSION
from ..utils.gregion import get_fixsize_regions, load_regions, output_regions

COMMAND = "convert"

def convert_main(argv):
    # parse command line options
    if len(argv) < 3:
        sys.stderr.write("Welcome to %s %s v%s!\n\n" % (APP, COMMAND, VERSION))
        sys.stderr.write("use -h or --help for help on argument.\n")
        sys.exit(1)

    parser = OptionParser(usage = "Usage: %s %s [options]" % (APP, COMMAND))
    parser.add_option("--input", "-i", dest="in_file", default=None,
        help=("Path to input region file."))
    parser.add_option("--inType", "-I", dest="in_type", default=None,
        help=("Input region type, one of bed|gff|tsv."))
    parser.add_option("--output", "-o", dest="out_file", default=None,
        help=("Path to output file; if not set, use stdout."))
    parser.add_option("--outType", "-O", dest="out_type", default="tsv",
        help=("Output region type, one of bed|tsv [default: %default]"))
    parser.add_option("--binsize", "-B", type="int", dest="bin_size", default=None,
        help=("Fixed size of bin in kb. it will be used when no input file."))
    parser.add_option("--hgver", "-H", type="int", dest="hg_ver", default=38,
        help=("Version of human genome, one of 19|38; set together with @p binsize [default: %default]"))

    # check options and args
    (options, args) = parser.parse_args(args = argv[2:])

    in_file, in_type = options.in_file, options.in_type
    out_file, out_type = options.out_file, options.out_type
    bin_size, hg_ver = options.bin_size, options.hg_ver

    if not in_file or not in_type:
        in_file = None
        in_type = None
    elif not os.path.isfile(in_file):
        sys.stderr.write("Error: input region file not exist: %s\n" % in_file)
        sys.exit(1)
    elif in_type.lower() not in ("bed", "gff", "tsv"):
        sys.stderr.write("Error: input region type should be one of bed|gff|tsv.\n")
        sys.exit(1)
    else:
        bin_size = None
        hg_ver = None
        in_type = in_type.lower()

    if in_file is None and (not bin_size or bin_size <= 0 or not hg_ver or hg_ver not in (19, 38)):
        sys.stderr.write("Error: either region file & type or a valid bin size & hg ver should be provided!\n")
        sys.exit(1)

    if not out_file:
        out_file = None

    if not out_type:
        sys.stderr.write("Error: out region type should be provided!\n")
        sys.exit(1)
    elif out_type not in ("bed", "tsv"):
        sys.stderr.write("Error: out region type should be one of bed|tsv.\n")
        sys.exit(1)

    # load regions
    reg_list = None
    if in_file is None:
        reg_list = get_fixsize_regions(bin_size, hg_ver)
    else:
        reg_list = load_regions(in_file, in_type)
    if not reg_list:
        sys.stderr.write("Error: empty region file or failed to parse regions.\n")
        sys.exit(1)

    # output regions
    output_regions(reg_list, out_file, out_type)