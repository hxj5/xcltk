#!/bin/bash
# baf_pre_phase.mouse.sh - call germline SNPs from mouse data.

# TODO: 
# - Test genotyping 
#   1) using common SNPs instead of performing denovo genotyping;
#   2) using normal cells only (or simply specify barcodes of normal cells in `-b`).
# - Test reference phasing locally instead of using online Sanger phasing.
# - Add reference version into the VCF header
# - Support calling germline SNPs for multiple bam files (-L)


function usage() {
    echo
    echo "Usage: $prog [options]"
    echo
    echo "Options:"
    echo "  -N, --name STR      Sample name."
    echo "  -s, --bam FILE      Path to bam file."
    echo "  -b, --barcode FILE  Path to barcode file. One barcode per line."
    echo "  -F, --phaseFA FILE  Path to mouse fasta file, e.g., mm10."
    echo "  -O, --outdir DIR    Path to output dir."
    echo "  -C, --celltag STR   Cell tag [${def_cell_tag}]"
    echo "  -u, --umi STR       UMI tag. Set to None to count reads [${def_umi}]"
    echo "  -D, --noDUP         If use, duplicate reads will be excluded."
    echo "  -p, --ncores INT    Number of cores [${def_ncores}]"
    echo "  -S, --smartseq      The input data is SMART-seq."
    echo "  -h, --help          Print this message and exit."
    echo
}


# global settings 
work_dir=`cd $(dirname $0); pwd`
prog=baf_pre_phase.mouse.sh

ucsc2ensembl=$work_dir/data/ucsc2ensembl.mouse.txt

# default settings
def_ncores=1
def_cell_tag=CB
def_umi=UB
use_dup=1
is_smart_seq=0


# check settings
if [ ! -e "$work_dir/utils.sh" ]; then
    echo "Error: utils file $work_dir/utils.sh does not exist!" >&2
    exit 1
fi
source $work_dir/utils.sh

assert_e  "$ucsc2ensembl"  "ucsc2ensembl file"


# parse args
if [ $# -lt 1 ]; then
    usage
    exit 1
fi

cmdline=`echo $0 $*`

ARGS=`getopt -o N:s:b:F:O:C:u:Dp:Sh --long name:,bam:,barcode:,phaseFA:,outdir:,celltag:,umi:,noDUP,ncores:,smartseq,help -n "" -- "$@"`
if [ $? -ne 0 ]; then
    log_err "Error: failed to parse command line. Terminating ..."
    exit 1
fi

eval set -- "$ARGS"
while true; do
    case "$1" in
        -N|--name) sid=$2; shift 2;;
        -s|--bam) bam=$2; shift 2;;
        -b|--barcode) barcode=$2; shift 2;;
        -F|--phaseFA) fa_phase=$2; shift 2;;
        -O|--outdir) out_dir=$2; shift 2;;
        -C|--celltag) cell_tag=$2; shift 2;;
        -u|--umi) umi=$2; shift 2;;
        -D|--noDUP) use_dup=0; shift;;
        -p|--ncores) ncores=$2; shift 2;;
        -S|--smartseq) is_smart_seq=1; shift;;
        -h|--help) usage; shift; exit 0;;
        --) shift; break;;
        *) log_err "Internal error!"; exit 1;;
    esac
done

log_msg "CMD: $cmdline"
set -x


# check cmdline args
assert_n  "$sid"  "Sample name"
assert_e  "$bam"  "BAM file"
assert_e  "$barcode"  "Barcode file"
assert_e  "$fa_phase"  "mouse fasta file"

assert_n  "$out_dir"  "Output dir"
if [ ! -e "$out_dir" ]; then
    mkdir -p $out_dir
fi
out_dir=`cd $out_dir; pwd`

if [ -z "$cell_tag" ]; then
    cell_tag=$def_cell_tag
fi

if [ -z "$umi" ]; then
    umi=$def_umi
fi

if [ $use_dup -eq 0 ]; then
    excl_flag=1796
else
    excl_flag=772
fi

if [ -z "$ncores" ]; then
    ncores=$def_ncores
fi


###### Core Part ######

res_dir=$out_dir/result
if [ ! -e "$res_dir" ]; then
    mkdir -p $res_dir
fi

# cell filtering
flt_bam=$res_dir/${sid}.cell_filter.bam
if [ $is_smart_seq -eq 0 ]; then
    log_msg "Filter cells. Only keep cells with valid barcodes."
    samtools view -h -b -D ${cell_tag}:${barcode} -@ $ncores $bam > $flt_bam
    samtools index $flt_bam
else
    log_msg "The input is SMART-seq data; Skip cell filtering."
    flt_bam=$bam
fi


# call germline SNPs
raw_vname=${sid}.raw.vcf.gz
chroms="`seq 1 19` X Y"
chroms=`echo $chroms | tr ' ' ',' | sed 's/,$//'`

log_msg "Call germline SNPs."

cellsnp-lite  -s $flt_bam  -O $res_dir/cellsnp            \
    --chrom $chroms  -p $ncores                            \
    --minMAF 0.1  --minCOUNT 20                            \
    --minLEN 30  --minMAPQ 20  --exclFLAG $excl_flag       \
    --cellTAG None  --UMItag $umi                          \
    --gzip  --genotype

raw_vpath=$res_dir/cellsnp/cellSNP.cells.vcf.gz


# Sanger Imputation Server fasta: chroms have no leading 'chr'
# bcftools annotate --rename-chrs is to removing the leading 'chr'
qc_vname=${raw_vname/.vcf/.het.qc.vcf}
qc_vpath=$res_dir/$qc_vname

log_msg "Keep heterozygous SNPs only."

bcftools view -Ou $raw_vpath |                          \
    bcftools view -Ou -i 'GT = "het"' |                  \
    bcftools view -Ou -i 'TYPE = "snp"' |                \
    bcftools annotate -Ou --rename-chrs $ucsc2ensembl |  \
    bcftools view -Oz -t $chroms > $qc_vpath


# filter by GQ
gq_bed=$res_dir/${qc_vname%.vcf.gz}.gq.bed
gq_vname=${qc_vname/.vcf/.gq.vcf}
gq_vpath=$res_dir/$gq_vname

log_msg "Filter SNPs with low GQ score."

bin_pl2gq=$res_dir/pl2gq.awk
generate_bin_pl2gq $bin_pl2gq

bcftools view -Ou $qc_vpath |                      \
    bcftools query -f '%CHROM\t%POS[\t%PL]\n' |    \
    $bin_pl2gq |                                   \
    awk '$NF > 20 { printf("%s\t%d\t%d\t%s\n", $1, $2 - 1, $2, $NF) }' > $gq_bed

bcftools view -Ou $qc_vpath |     \
    bcftools view -Oz -T $gq_bed > $gq_vpath

flt_vname=$gq_vname
flt_vpath=$gq_vpath


# bcftools fixref check
log_msg "(Pre-xcltk-fixref) bcftools fixref check."
bcftools +fixref $flt_vpath -- -f $fa_phase


# xcltk fixref
fix_vname=${flt_vname/.vcf/.fixref.sort.vcf}
fix_vpath=$res_dir/$fix_vname

log_msg "xcltk fixref."

xcltk fixref  -i $flt_vpath  -r $fa_phase  -v |   \
    bcftools sort -Oz > $fix_vpath


# bcftools fixref check
log_msg "(Post-xcltk-fixref) bcftools fixref check."
bcftools +fixref $fix_vpath -- -f $fa_phase


# move final result
mv $fix_vpath $out_dir


###### END ######
log_msg "All Done!"
log_msg "End"

