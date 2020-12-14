#!/bin/bash
#install all softwares needed for xcltk pipeline

if [ $# -lt 2 ]; then
    echo "Usage: $0 <out_dir> <conda env>" >&2
    exit 1
fi

out_dir=$1
env=$2

work_dir=`cd $(dirname $0); pwd`
app_lst=$work_dir/conda.packages.all.lst
if [ ! -f "$app_lst" ]; then
    echo "Conda package list file 'conda.packages.all.lst' needed!" >&2
    exit 3
fi


