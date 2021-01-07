#!/bin/bash
#this script is aimed to build a pipeline for post-imputation.

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

# other scripts
bin_pre_pileup=$work_dir/baf_pre_pileup.sh
bin_pileup=$work_dir/baf_pileup.sh

###### Init running ######
# import utils
source $utils        # import eval_cmd, load_cfg, log_msg, log_err

function usage() {
    echo
    echo "Usage: $1 [options]"
    echo
    echo "Options:"
    echo "  -N, --name STR       Sample name"
    echo "  -S, --seq STR        Seq type: dna|rna|atac"
    echo "  -b, --barcode FILE   Path to barcode file, one barcode per line"
    echo "  -s, --bam FILE       Path to bam file for droplet-based dataset"
    echo "  -L, --bamlist File   Path to bam list file for well-based dataset"
    echo "  -u, --umi STR        UMI tag if available"
    echo "  -f, --fasta FILE     Path to fasta file"
    echo "  -g, --hg INT         Version of fasta, 19 or 38"
    echo "  -P, --phase STR      Phase type: phase|impute"
    echo "  -v, --vcf FILE       Path to phased vcf"
    echo "  -B, --blocks FILE    Region of feature blocks in bed|tsv|gff"
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

ARGS=`getopt -o N:S:P:s:L:u:f:g:b:v:B:p:O:c:h --long name:,seq:,phase:,bam:,bamlist:,umi:,fasta:,hg:,barcode:,vcf:,blocks:,ncores:,outdir:,config:,help -n "" -- "$@"`
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
        -s|--bam) bam=$2; shift 2;;
        -L|--bamlist) bam_list=$2; shift 2;;
        -u|--umi) umi=$2; shift 2;;
        -f|--fasta) fasta=$2; shift 2;;
        -g|--hg) hg=$2; shift 2;;
        -b|--barcode) barcode=$2; shift 2;;
        -v|--vcf) vcf=$2; shift 2;;
        -B|--blocks) blocks=$2; shift 2;;
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

if [ -n "$bam" ]; then            # droplet-based dataset
    bam_opt="-s $bam -b $barcode"
elif [ -n "$bam_list" ]; then     # well-based dataset
    bam_opt="-S $bam_list"
else
    log_err "Error: either --bam or --bamlist should be specified!"
    exit 1
fi

if [ -z "$umi" ]; then
    umi_opt=""
else
    umi_opt="-u $umi"
fi

sid=${smp_name}

###### Core Part ######
aim="run pre-pileup pipeline"
cmd="$bin_pre_pileup -N $smp_name -S $seq_type -P $phase_type -f $fasta \\
       -g $hg -v $vcf -O $out_dir -c $cfg"
eval_cmd "$cmd" "$aim"

preplp_vpath=`ls $out_dir/*.uniq.sort.vcf.gz`

aim="run pileup"
cmd="$bin_pileup -S $seq_type $bam_opt $umi_opt -v $preplp_vpath \\
       -p $ncores -O $out_dir -c $cfg"
eval_cmd "$cmd" "$aim"

xcsp_dir=$out_dir/xcltk-pileup
gt_vpath=`ls $out_dir/*.gt.vcf.gz`

aim="phase SNPs into haplotype blocks of even size"
phs_even_dir=$out_dir/phase-snp-even
mkdir -p $phs_even_dir &> /dev/null
size=50    # kb
cmd="$bin_xcltk convert -B $size -H $hg -o $phs_even_dir/blocks.${size}kb.tsv && \\
     $bin_xcltk phase_snp --sid $sid --snpAD $xcsp_dir/cellSNP.tag.AD.mtx \\
       --snpDP $xcsp_dir/cellSNP.tag.DP.mtx --phase $gt_vpath \\
       --region $phs_even_dir/blocks.${size}kb.tsv --outdir $phs_even_dir"
eval_cmd "$cmd" "$aim"

cp $out_dir/cellsnp-lite/cellSNP.samples.tsv $phs_even_dir

aim="phase SNPs into haplotype blocks of features"
phs_fet_dir=$out_dir/phase-snp-feature
mkdir -p $phs_fet_dir &> /dev/null
cmd="$bin_xcltk phase_snp --sid $sid --snpAD $xcsp_dir/cellSNP.tag.AD.mtx \\
       --snpDP $xcsp_dir/cellSNP.tag.DP.mtx \\
       --phase $gt_vpath --region $blocks --outdir $phs_fet_dir"
eval_cmd "$cmd" "$aim"

cp $out_dir/cellsnp-lite/cellSNP.samples.tsv $phs_fet_dir
cp $blocks $phs_fet_dir

###### END ######
log_msg "All Done!"
log_msg "End"

