#!/bin/bash
#install all softwares needed by xcltk preprocessing pipeline

if [ $# -lt 2 ]; then
    echo >&2
    echo "Note: This script is aimed to install all softwares needed by the" >&2
    echo "      preprocessing pipeline to certain conda environment." >&2
    echo >&2
    echo "Usage: $0 <out dir> <conda env>" >&2
    echo >&2
    exit 1
fi

out_dir=$1
env=$2

work_dir=`cd $(dirname $0); pwd`

# 
pkg_lst=$work_dir/preprocess.conda.packages.lst
if [ ! -f "$pkg_lst" ]; then
    echo "Error: conda package list file 'preprocess.conda.packages.lst' needed!" >&2
    exit 3
fi

