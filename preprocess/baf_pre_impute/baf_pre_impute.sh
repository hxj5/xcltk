#!/bin/bash
#this script is aimed to build a pipeline for steps from calling germline SNPs to pre-imputation.

###### Global settings ######
work_dir=`cd $(dirname $0); pwd`
prog_path=$0
prog_name=`basename $0`

# default config file 
cfg=$work_dir/${prog_name%.sh}.cfg

# utils file
utils=$work_dir/../utils.sh

# ucsc2ensembl file
ucsc2ensembl=$work_dir/../annotate/ucsc2ensembl.txt

# some scripts
bin_gl2gq=$work_dir/gl2gq.awk
bin_pl2gq=$work_dir/pl2gq.awk

# liftover files
bin_py_liftover=$work_dir/../liftover/liftOver_vcf.py
chain_hg38to19=$work_dir/../liftover/hg38ToHg19.over.chain.gz

###### Init running ######
# import utils
source $utils        # import eval_cmd, load_cfg, log_msg, log_err

###### Command Line Parsing ######
function usage() {
    echo
    echo "Usage: $1 [options]"
    echo
    echo "Options:"
    echo "  -N, --name STR      Sample name"
    echo "  -s, --bam FILE      Path to bam file"
    echo "  -f, --fasta FILE    Path to fasta file"
    echo "  -g, --hg INT        Version of fasta, 19 or 38"
    echo "  -O, --outdir DIR    Path to output dir"
    echo "  -c, --config FILE   Path to config file. If not set, use the"
    echo "                      default config baf_pre_impute.cfg"
    echo "  -h, --help          This message"
    echo
    echo "Note:"
    echo "  It's recommended to make a copy of the config file baf_pre_impute.cfg"
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

ARGS=`getopt -o N:s:f:g:O:c:h --long name:,bam:,fasta:,hg:,outdir:,config:,help -n "" -- "$@"`
if [ $? -ne 0 ]; then
    log_err "Error: failed to parse command line. Terminating..."
    exit 1
fi

eval set -- "$ARGS"
while true; do
    case "$1" in
        -N|--name) smp_name=$2; shift 2;;
        -s|--bam) bam=$2; shift 2;;
        -f|--fasta) fasta=$2; shift 2;;
        -g|--hg) hg=$2; shift 2;;
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

sid=$smp_name       # Sample ID

mkdir -p $out_dir &> /dev/null
out_dir=`cd $out_dir; pwd`

if [ -n "$hg" ]; then
    if [ $hg -ne 19 ] && [ $hg -ne 38 ]; then
        log_err "Error: hg version should be 19 or 38"
        exit 1
    fi
else
    log_err "Error: hg version empty!"
    exit 1
fi

###### Core Part ######
aim="call germline SNPs"
raw_vname=${sid}.hg${hg}.raw.vcf.gz
raw_vpath=$out_dir/${sid}.hg${hg}.raw.vcf.gz
cmd="$bin_freebayes -C 2 -F 0.1 -m 20 --min-coverage 20 -f $fasta $bam | \\
     $bin_bgzip -c > $raw_vpath"
eval_cmd "$cmd" "$aim"

# Sanger Imputation Server fasta: chroms have no leading 'chr'
# bcftools annotate --rename-chrs is to removing the leading 'chr'
aim="QC: \\
      + filter low QUAL and DP;                                          \\
      + filter records that are not of SNP type;                         \\
      + rename chroms, remove the leading 'chr' from the name of chroms; \\
      + filter records not in target chroms (default chr1-22, X, Y);"
qc_vname=${raw_vname/.vcf/.qc.vcf}
qc_vpath=$out_dir/$qc_vname
target_chroms="`seq 1 22` X Y"
tgt_chroms=`echo $target_chroms | tr ' ' ',' | sed 's/,$//'`
cmd="$bin_bcftools view -Ou $raw_vpath |                         \\
     $bin_bcftools view -Ou -i 'QUAL > 20 && INFO/DP > 0' |       \\
     $bin_bcftools view -Ou -i 'TYPE = \"snp\"' |                  \\
     $bin_bcftools annotate -Ou --rename-chrs $ucsc2ensembl |     \\
     $bin_bcftools view -Oz -t $tgt_chroms                        \\
      > $qc_vpath"
eval_cmd "$cmd" "$aim"

aim="filter by GQ"
gq_bed=${qc_vname%.vcf.gz}.gq.bed
gq_vname=${qc_vname/.vcf/.gq.vcf}
gq_vpath=$out_dir/$gq_vname
cmd="$bin_bcftools view -Ou $qc_vpath |                          \\
     $bin_bcftools query -f '%CHROM\t%POS[\t%GL]\n' |             \\
     $bin_gl2gq |                                                \\
     awk '\$NF > 20 { printf(\"%s\t%d\t%d\t%s\n\", \$1, \$2 - 1, \$2, \$NF) }' \\
       > $gq_bed &&                                                            \\
     $bin_bcftools view -Ou $qc_vpath |                                     \\
     $bin_bcftools view -Oz -T $gq_bed                                      \\
       > $gq_vpath"
eval_cmd "$cmd" "$aim"

flt_vname=$gq_vname
flt_vpath=$gq_vpath

aim="convert hg38 to hg19"
lift_vname=${flt_vname/.vcf/.hg19.vcf}
lift_vpath=$out_dir/$lift_vname
if [ $hg -eq 19 ]; then
    log_msg "$aim"
    log_msg "Already hg19, skip liftover"
    lift_vname=$flt_vname
    lift_vpath=$flt_vpath
else
    cmd="$bin_python $bin_py_liftover -c $chain_hg38to19 -i $flt_vpath \\
           -o ${lift_vpath}.tmp -P $bin_liftover &&                    \\
         $bin_bcftools view -i 'POS > 0' -Oz ${lift_vpath}.tmp         \\
           > ${lift_vpath} &&                                          \\
         rm ${lift_vpath}.tmp"
    eval_cmd "$cmd" "$aim"
fi

aim="bcftools fixref checking"
cmd="$bin_bcftools +fixref $lift_vpath -- -f $fa_impute"
eval_cmd "$cmd" "$aim"

aim="xcltk fixref"
fix_vname=${lift_vname/.vcf/.fixref.sort.vcf}
fix_vpath=$out_dir/$fix_vname
cmd="$bin_bcftools query -f '%CHROM:%POS-%POS\n' $lift_vpath          \\
       > ${lift_vpath}.region.lst &&                                  \\
     $bin_samtools faidx -r ${lift_vpath}.region.lst $fa_impute |     \\
     $bin_bgzip -c > ${lift_vpath}.fa.gz &&                           \\
     $bin_xcltk fixref -i $lift_vpath -r ${lift_vpath}.fa.gz |        \\
     $bin_bcftools sort -Oz > $fix_vpath &&                           \\
     rm ${lift_vpath}.region.lst"
eval_cmd "$cmd" "$aim"

aim="bcftools fixref checking"
cmd="$bin_bcftools +fixref $fix_vpath -- -f $fa_impute"
eval_cmd "$cmd" "$aim"

###### END ######
log_msg "All Done!"
log_msg "End"

