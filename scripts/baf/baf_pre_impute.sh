#!/bin/bash
#this script is aimed to build a pipeline for steps from calling germline SNPs to pre-imputation.

bin_bcftools=bcftools
bin_cellsnp_lite=cellsnp-lite
bin_freebayes=freebayes
bin_python=python
bin_samtools=samtools

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


