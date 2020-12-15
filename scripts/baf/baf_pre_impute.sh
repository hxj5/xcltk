#!/bin/bash
#this script is aimed to build a pipeline for steps from calling germline SNPs to pre-imputation.

###### Global settings ######
work_dir=`cd $(dirname $0); pwd`

# config file for binary paths
bin_cfg=$work_dir/baf_pre_impute.bin.cfg

# utils file
utils=$work_dir/../utils.sh

function print_usage() {
    echo
    echo "Usage: $0 [options]"
    echo
    echo "Options:"
    echo "  -N, --name STR    Sample name."
    echo "  -s, --bam FILE    Path to bam file."
    echo "  -O, --outdir DIR  Path to output dir."
    echo "  --rootdir DIR     Path to root dir of this project."
    echo "  --remove-tmp      If set, tmp files will be removed."
    echo "  -h, --help        This message."
    echo
}


