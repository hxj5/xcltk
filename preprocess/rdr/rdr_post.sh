#!/bin/bash

if [ $# -lt 1 ]; then
    echo "Usage $0 <rdr_dir>" >&2
    exit 1
fi

rdr_dir=$1

# rename matrix to rdr.mtx
if [ -e "$rdr_dir/matrix.mtx" ]; then
    mv $rdr_dir/matrix.mtx $rdr_dir/rdr.mtx
fi
