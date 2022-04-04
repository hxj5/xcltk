#!/bin/bash
#this script is aimed to build a pipeline for pileup and post_pileup.

# TODO: add reference version into VCF header

###### Global settings ######
work_dir=`cd $(dirname $0); pwd`
prog_path=$0
prog_name=`basename $0`

# default config file 
cfg_fn=${prog_name%.sh}.cfg
cfg=$work_dir/$cfg_fn

# utils file
utils=$work_dir/../utils/utils.sh

###### Init running ######
# import utils
source $utils        # import eval_cmd, load_cfg, log_msg, log_err

function usage() {
    echo
    echo "Usage: $1 [options]"
    echo
    echo "Options:"
    echo "  -S, --seq STR        Seq type: dna|rna|atac"
    echo "  -b, --barcode FILE   Path to barcode file, one barcode per line"
    echo "  -s, --bam FILE       Path to bam file for droplet-based dataset"
    echo "  -L, --bamlist File   Path to bam list file for well-based dataset"
    echo "  -u, --umi STR        UMI tag if available"
    echo "  -B, --blocks FILE    Region of feature blocks in TSV format"
    echo "  -v, --vcf FILE       Path to phased vcf"
    echo "  -g, --hg INT         Version of fasta, 19 or 38"
    echo "  -p, --ncores INT     Number of cores"
    echo "  -O, --outdir DIR     Path to output dir"
    echo "  -c, --config FILE    Path to config file. If not set, use the"
    echo "                       default config $cfg_fn"
    echo "  -h, --help           This message"
    echo
    echo "Note:"
    echo "  It's recommended to make a copy of the config file $cfg_fn"
    echo "  and then modify the new copy instead of modifying the original file"
    echo "  in place."
    echo
}

cmdline=`echo $0 $*`
log_msg "Start ..."
log_msg "Command Line: $cmdline"

# parse args
log_msg "Parse cmdline ..."
if [ $# -lt 1 ]; then
    usage $prog_path
    exit 1
fi

ARGS=`getopt -o S:s:b:L:u:B:v:g:p:O:c:h --long seq:,bam:,barcode:,bamlist:,umi:,blocks:,vcf:,hg:,ncores:,outdir:,config:,help -n "" -- "$@"`
if [ $? -ne 0 ]; then
    echo "Error: failed to parse command line args. Terminating..." >&2
    exit 1
fi

eval set -- "$ARGS"
while true; do
    case "$1" in
        -S|--seq) seq_type=$2; shift 2;;
        -s|--bam) bam=$2; shift 2;;
        -b|--barcode) barcode=$2; shift 2;;
        -L|--bamlist) bam_list=$2; shift 2;;
        -u|--umi) umi=$2; shift 2;;
        -B|--blocks) blocks_fet=$2; shift 2;;
        -v|--vcf) vcf=$2; shift 2;;
        -g|--hg) hg=$2; shift 2;;
        -p|--ncores) ncores=$2; shift 2;;
        -O|--outdir) out_dir=$2; shift 2;;
        -c|--config) cfg=$2; shift 2;;
        -h|--help) usage $prog_path; shift; exit 0;;
        --) shift; break;;
        *) log_err "Internal error!"; exit 1;;
    esac
done

# check cmdline args
# TODO: check each arg, is null or path exists?
log_msg "Load config ..."
load_cfg $cfg

mkdir -p $out_dir &> /dev/null
out_dir=`cd $out_dir; pwd`

if [ -n "$bam" ] && [ -n "$bam_list" ]; then
    log_err "Error: --bam and --bamlist should not be specified at the same time!"
    exit 1
fi

cell_tag=CB
if [ -n "$bam" ]; then            # droplet-based dataset
    bam_opt="-s $bam -b $barcode"
elif [ -n "$bam_list" ]; then     # well-based dataset
    bam_opt="-S $bam_list"
    cell_tag=None
else
    log_err "Error: either --bam or --bamlist should be specified!"
    exit 1
fi

if [ -z "$umi" ]; then 
    umi=None
    excl_flag=1796
else
    excl_flag=772
fi

csp_in_vpath=$vcf
csp_in_vname=`basename $vcf`

###### Core Part ######
aim="cellsnp-lite pileup"
csp_dir=$out_dir/cellsnp-lite
if [ "$seq_type" == "dna" ] || [ "$seq_type" == "atac" ]; then
    cmd="$bin_cellsnp $bam_opt -O $csp_dir -R $csp_in_vpath --minCOUNT 1 --minMAF 0 \\
      --minLEN 30 --minMAPQ 20 --inclFLAG 0 --exclFLAG $excl_flag --UMItag $umi -p $ncores  \\
      --cellTAG $cell_tag --genotype --gzip"
elif [ "$seq_type" == "rna" ]; then
    cmd="$bin_cellsnp $bam_opt -O $csp_dir -R $csp_in_vpath --minCOUNT 1 --minMAF 0 \\
      --minLEN 30 --minMAPQ 20 --inclFLAG 0 --exclFLAG $excl_flag --UMItag $umi -p $ncores  \\
      --cellTAG $cell_tag --genotype --gzip"
else  # unknown
    log_msg "Warning: unknown seq type $seq_type, use the dna pileup method"
    cmd="$bin_cellsnp $bam_opt -O $csp_dir -R $csp_in_vpath --minCOUNT 1 --minMAF 0  \\
      --minLEN 30 --minMAPQ 20 --inclFLAG 0 --exclFLAG $excl_flag --UMItag $umi -p $ncores   \\
      --cellTAG $cell_tag --genotype --gzip"
fi
eval_cmd "$cmd" "$aim"

csp_vpath=$csp_dir/cellSNP.base.vcf.gz

aim="merge pileup vcf and phase GT vcf"
gt_vname=${csp_in_vname/.vcf/.gt.vcf}
gt_vpath=$out_dir/$gt_vname
cmd="zcat $csp_in_vpath | sed 's/^chr//' | $bin_bgzip -c > ${csp_in_vpath}.tmp &&
     zcat $csp_vpath | sed 's/^chr//' | $bin_bgzip -c > ${csp_vpath}.tmp &&
     $bin_bcftools view -Oz -T ${csp_vpath}.tmp ${csp_in_vpath}.tmp > $gt_vpath &&
     rm ${csp_in_vpath}.tmp && rm ${csp_vpath}.tmp"
eval_cmd "$cmd" "$aim"

aim="create file containing blocks of even size"
phs_even_dir=$out_dir/phase-snp-even
mkdir -p $phs_even_dir &> /dev/null
size=50    # kb
blocks_even=$phs_even_dir/blocks.${size}kb.tsv
cmd="$bin_xcltk convert -B $size -H $hg -o ${blocks_even}.tmp && 
     cat ${blocks_even}.tmp | awk '{printf(\"%s\t%s:%s-%s\", \$0, \$1, \$2, \$3)}' > ${blocks_even} &&
     rm ${blocks_even}.tmp"
eval_cmd "$cmd" "$aim"

aim="phase SNPs into haplotype blocks of even size"
if [ "$seq_type" == "dna" ] || [ "$seq_type" == "atac" ]; then
    cmd="$bin_xcltk pileup $bam_opt -O $phs_even_dir -R $blocks_even -P $gt_vpath --minCOUNT 1 --minMAF 0 \\
      --minLEN 30 --minMAPQ 20 --inclFLAG 0 --exclFLAG $excl_flag --UMItag $umi -p $ncores \\
      --cellTAG $cell_tag"
elif [ "$seq_type" == "rna" ]; then
    cmd="$bin_xcltk pileup $bam_opt -O $phs_even_dir -R $blocks_even -P $gt_vpath --minCOUNT 1 --minMAF 0 \\
      --minLEN 30 --minMAPQ 20 --inclFLAG 0 --exclFLAG $excl_flag --UMItag $umi -p $ncores \\
      --cellTAG $cell_tag"
else  # unknown
    log_msg "Warning: unknown seq type $seq_type, use the dna pileup method"
    cmd="$bin_xcltk pileup $bam_opt -O $phs_even_dir -R $blocks_even -P $gt_vpath --minCOUNT 1 --minMAF 0 \\
      --minLEN 30 --minMAPQ 20 --inclFLAG 0 --exclFLAG $excl_flag --UMItag $umi -p $ncores \\
      --cellTAG $cell_tag"
fi
eval_cmd "$cmd" "$aim"

aim="phase SNPs into haplotype blocks of features"
phs_fet_dir=$out_dir/phase-snp-feature
mkdir -p $phs_fet_dir &> /dev/null
if [ "$seq_type" == "dna" ] || [ "$seq_type" == "atac" ]; then
    cmd="$bin_xcltk pileup $bam_opt -O $phs_fet_dir -R $blocks_fet -P $gt_vpath --minCOUNT 1 --minMAF 0 \\
      --minLEN 30 --minMAPQ 20 --inclFLAG 0 --exclFLAG $excl_flag --UMItag $umi -p $ncores \\
      --cellTAG $cell_tag"
elif [ "$seq_type" == "rna" ]; then
    cmd="$bin_xcltk pileup $bam_opt -O $phs_fet_dir -R $blocks_fet -P $gt_vpath --minCOUNT 1 --minMAF 0 \\
      --minLEN 30 --minMAPQ 20 --inclFLAG 0 --exclFLAG $excl_flag --UMItag $umi -p $ncores \\
      --cellTAG $cell_tag"
else  # unknown
    log_msg "Warning: unknown seq type $seq_type, use the dna pileup method"
    cmd="$bin_xcltk pileup $bam_opt -O $phs_fet_dir -R $blocks_fet -P $gt_vpath --minCOUNT 1 --minMAF 0 \\
      --minLEN 30 --minMAPQ 20 --inclFLAG 0 --exclFLAG $excl_flag --UMItag $umi -p $ncores \\
      --cellTAG $cell_tag"
fi
eval_cmd "$cmd" "$aim"

###### END ######
log_msg "All Done!"
log_msg "End"
