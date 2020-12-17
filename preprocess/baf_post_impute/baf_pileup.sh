#!/bin/bash
#this script is aimed to build a pipeline for steps from post-imputation to pileup.

###### Global settings ######
work_dir=`cd $(dirname $0); pwd`
prog_path=$0
prog_name=`basename $0`

# default config file 
cfg=$work_dir/${prog_name%.sh}.cfg

# utils file
utils=$work_dir/../utils.sh

# ensembl2ucsc file
=$work_dir/impute/ucsc2ensembl.txt

function usage() {
    echo
    echo "Usage: $1 [options]"
    echo
    echo "Options:"
    echo "  -N, --name STR     Sample name."
    echo "  -S, --seq STR      Seq type: dna|rna|atac."
    echo "  -P, --phase STR    Phase type: phase|impute."
    echo "  -s, --bam FILE     Path to bam file."
    echo "  -b, --barcode FILE Path to barcode file."
    echo "  -v, --vcf FILE     Path to phased vcf."
    echo "  -p, --ncores INT   Number of cores."
    echo "  -O, --outdir DIR   Path to output dir."
    echo "  -h, --help         This message."
    echo
}

### parse command line args
if [ $# -lt 1 ]; then
    print_usage $script_name
    exit 1
fi

# parse args
cmdline=`echo $0 $*`
is_remove_tmp=0
ARGS=`getopt -o N:S:P:s:b:v:p:O:h --long name:,seq:,phase:,bam:,barcode:,vcf:,ncores:,rootdir:,outdir:,remove-tmp,help -n "" -- "$@"`
if [ $? -ne 0 ]; then
    echo "Error: failed to parse command line args. Terminating..." >&2
    exit 1
fi
eval set -- "$ARGS"
while true; do
    case "$1" in
        -N|--name) smp_name=$2; shift 2;;
        -S|--seq) seq_type=$2; shift 2;;
        -P|--phase) phase_type=$2; shift 2;;
        -s|--bam) bam_file=$2; shift 2;;
        -b|--barcode) barcode_file=$2; shift 2;;
        -v|--vcf) vcf_file=$2; shift 2;;
        -p|--ncores) ncores=$2; shift 2;;
        -O|--outdir) out_dir=$2; shift 2;;
        --rootdir) root_dir=$2; shift 2;;
        --remove-tmp) is_remove_tmp=1; shift;;
        -h|--help) print_usage $script_name; shift; exit 0;;
        --) shift; break;;
        *) echo "Internal error!" >&2; exit 1;;
    esac
done

# check cmdline args
if [ -z $root_dir ] || [ ! -d $root_dir ]; then
    echo "Error: root_dir invalid!" >&2
    exit 1
fi
