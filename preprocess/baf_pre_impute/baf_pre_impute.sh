#!/bin/bash
#this script is aimed to build a pipeline for steps from calling germline SNPs to pre-imputation.

# TODO: 
# - add reference version into the VCF header
# - support calling germline SNPs for multiple bam files (-L)

###### Global settings ######
work_dir=`cd $(dirname $0); pwd`
prog_path=$0
prog_name=`basename $0`

# default config file 
cfg_fn=${prog_name%.sh}.cfg
cfg=$work_dir/$cfg_fn

# utils file
utils=$work_dir/../utils/utils.sh

# ucsc2ensembl file
ucsc2ensembl=$work_dir/../data/annotate/ucsc2ensembl.txt

# some scripts
bin_gl2gq=$work_dir/gl2gq.awk
bin_pl2gq=$work_dir/pl2gq.awk

# liftover files
bin_py_liftover=$work_dir/../utils/liftOver_vcf.py
chain_hg38to19=$work_dir/../data/liftover/hg38ToHg19.over.chain.gz

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
    echo "  -s, --bam FILE      Path to bam file for droplet-based dataset"
#   echo "  -L, --bamlist File  Path to bam list file for well-based dataset"
    echo "  -f, --fasta FILE    Path to fasta file"
    echo "  -g, --hg INT        Version of fasta, 19 or 38"
    echo "  -C, --call STR      Germline SNP calling app, cellsnp-lite or"
    echo "                      freebayes [cellsnp-lite]"
    echo "  -O, --outdir DIR    Path to output dir"
    echo "  -p, --ncores INT    Number of cores"
    echo "  -c, --config FILE   Path to config file. If not set, use the"
    echo "                      default config $cfg_fn"
    echo "  -h, --help          This message"
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

ARGS=`getopt -o N:s:L:f:g:C:O:p:c:h --long name:,bam:,bamList:,fasta:,hg:,call:,outdir:,ncores:,config:,help -n "" -- "$@"`
if [ $? -ne 0 ]; then
    log_err "Error: failed to parse command line. Terminating..."
    exit 1
fi

eval set -- "$ARGS"
while true; do
    case "$1" in
        -N|--name) smp_name=$2; shift 2;;
        -s|--bam) bam=$2; shift 2;;
        -L|--bamlist) bam_list=$2; shift 2;;
        -f|--fasta) fasta=$2; shift 2;;
        -g|--hg) hg=$2; shift 2;;
        -C|--call) app_call=$2; shift 2;;
        -O|--outdir) out_dir=$2; shift 2;;
        -p|--ncores) ncores=$2; shift 2;;
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

if [ -n "$bam" ] && [ -n "$bam_list" ]; then
    log_err "Error: --bam and --bamlist should not be specified at the same time!"
    exit 1
fi

if [ -n "$bam" ]; then            # droplet-based dataset
    bam_opt="-s $bam"
elif [ -n "$bam_list" ]; then     # well-based dataset
    bam_opt="-S $bam_list"
    if [ "$app_call" == "freebayes" ]; then
        log_err "Error: should not specify --bamList for freebayes!"
        exit 1
    fi
else
    log_err "Error: either --bam or --bamlist should be specified!"
    exit 1
fi

if [ -n "$hg" ]; then
    if [ $hg -ne 19 ] && [ $hg -ne 38 ]; then
        log_err "Error: hg version should be 19 or 38"
        exit 1
    fi
else
    log_err "Error: hg version empty!"
    exit 1
fi

if [ -z "$ncores" ]; then 
  ncores=1
fi

###### Core Part ######
aim="call germline SNPs"
raw_vname=${sid}.hg${hg}.raw.vcf.gz
raw_vpath=$out_dir/$raw_vname
target_chroms="`seq 1 22` X Y"
tgt_chroms=`echo $target_chroms | tr ' ' ',' | sed 's/,$//'`
if [ "$app_call" == "freebayes" ]; then
    cmd="$bin_freebayes -C 2 -F 0.1 -m 20 --min-coverage 20 -f $fasta $bam | 
         $bin_bgzip -c > $raw_vpath"
else
    cmd="$bin_cellsnp $bam_opt -O $out_dir/cellsnp_pre -p $ncores --minMAF 0.1 \\
         --minCOUNT 20 --minLEN 30 --minMAPQ 20 --exclFLAG 1796 --cellTAG None \\
         --UMItag None --chrom $tgt_chroms --gzip --genotype"
    raw_vpath=$out_dir/cellsnp_pre/cellSNP.cells.vcf.gz
fi
eval_cmd "$cmd" "$aim"

# Sanger Imputation Server fasta: chroms have no leading 'chr'
# bcftools annotate --rename-chrs is to removing the leading 'chr'
aim="QC: 
  + filter low QUAL and DP;                                          
  + filter records that are not of SNP type;                         
  + filter strlen(REF) != 1 || N_ALT != 1;
  + rename chroms, remove the leading 'chr' from the name of chroms; 
  + filter records not in target chroms (default chr1-22, X, Y);"
qc_vname=${raw_vname/.vcf/.qc.vcf}
qc_vpath=$out_dir/$qc_vname
if [ "$app_call" == "freebayes" ]; then
    cmd="$bin_bcftools view -Ou $raw_vpath | 
      $bin_bcftools view -Ou -i 'QUAL > 20 && INFO/DP > 0' |     
      $bin_bcftools view -Ou -i 'TYPE = \"snp\"' |        
      $bin_bcftools view -Ou -i 'STRLEN(REF) == 1 && N_ALT == 1' |
      $bin_bcftools annotate -Ou --rename-chrs $ucsc2ensembl | 
      $bin_bcftools view -Oz -t $tgt_chroms > $qc_vpath"
else
    cmd="$bin_bcftools view -Ou $raw_vpath | 
      $bin_bcftools view -Ou -i 'TYPE = \"snp\"' |        
      $bin_bcftools annotate -Ou --rename-chrs $ucsc2ensembl | 
      $bin_bcftools view -Oz -t $tgt_chroms > $qc_vpath"
fi
eval_cmd "$cmd" "$aim"

aim="filter by GQ"
gq_bed=$out_dir/${qc_vname%.vcf.gz}.gq.bed
gq_vname=${qc_vname/.vcf/.gq.vcf}
gq_vpath=$out_dir/$gq_vname
if [ "$app_call" == "freebayes" ]; then
    cmd="$bin_bcftools view -Ou $qc_vpath |                          
      $bin_bcftools query -f '%CHROM\t%POS[\t%GL]\n' |             
      $bin_gl2gq |                                                
      awk '\$NF > 20 { printf(\"%s\t%d\t%d\t%s\n\", \$1, \$2 - 1, \$2, \$NF) }' > $gq_bed &&                                                         
      $bin_bcftools view -Ou $qc_vpath |                                 
      $bin_bcftools view -Oz -T $gq_bed > $gq_vpath"
else
    cmd="$bin_bcftools view -Ou $qc_vpath |                          
      $bin_bcftools query -f '%CHROM\t%POS[\t%PL]\n' |             
      $bin_pl2gq |                                                
      awk '\$NF > 20 { printf(\"%s\t%d\t%d\t%s\n\", \$1, \$2 - 1, \$2, \$NF) }' > $gq_bed &&                                                         
      $bin_bcftools view -Ou $qc_vpath |                                 
      $bin_bcftools view -Oz -T $gq_bed > $gq_vpath"
fi
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
      -o ${lift_vpath/.vcf/.tmp.vcf} -P $bin_liftover &&  
      $bin_bcftools view -i 'POS > 0' -Oz ${lift_vpath/.vcf/.tmp.vcf} > ${lift_vpath} &&                                      
      rm ${lift_vpath/.vcf/.tmp.vcf}"
    eval_cmd "$cmd" "$aim"
fi

aim="bcftools fixref checking"
cmd="$bin_bcftools +fixref $lift_vpath -- -f $fa_impute"
eval_cmd "$cmd" "$aim"

aim="xcltk fixref"
fix_vname=${lift_vname/.vcf/.fixref.sort.vcf}
fix_vpath=$out_dir/$fix_vname
tmp_prefix=${lift_vpath%.vcf.gz}
cmd="$bin_bcftools query -f '%CHROM:%POS-%POS\n' $lift_vpath > ${tmp_prefix}.region.lst &&                                
  $bin_samtools faidx -r ${tmp_prefix}.region.lst $fa_impute |  
  $bin_bgzip -c > ${tmp_prefix}.fa.gz &&                       
  $bin_xcltk fixref -i $lift_vpath -r ${tmp_prefix}.fa.gz |   
  $bin_bcftools sort -Oz > $fix_vpath && 
  rm ${tmp_prefix}.region.lst"
eval_cmd "$cmd" "$aim"

aim="bcftools fixref checking"
cmd="$bin_bcftools +fixref $fix_vpath -- -f $fa_impute"
eval_cmd "$cmd" "$aim"

###### END ######
log_msg "All Done!"
log_msg "End"

