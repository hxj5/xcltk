#!/bin/bash
#this script is aimed to build a pipeline for steps from post-imputation to pileup.

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

# ensembl2ucsc file
ensembl2ucsc=$work_dir/../data/annotate/ensembl2ucsc.txt

# liftover files
bin_py_liftover=$work_dir/../utils/liftOver_vcf.py
chain_hg19to38=$work_dir/../data/liftover/hg19ToHg38.over.chain.gz

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
    echo "  -b, --barcode FILE   Path to barcode file"
    echo "  -s, --bam FILE       Path to bam file, one barcode per line"
    echo "  -f, --fasta FILE     Path to fasta file"
    echo "  -g, --hg INT         Version of fasta, 19 or 38"
    echo "  -P, --phase STR      Phase type: phase|impute"
    echo "  -v, --vcf FILE       Path to phased vcf"
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

ARGS=`getopt -o N:S:P:s:f:g:b:v:p:O:c:h --long name:,seq:,phase:,bam:,fasta:,hg:,barcode:,vcf:,ncores:,outdir:,config:,help -n "" -- "$@"`
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
        -f|--fasta) fasta=$2; shift 2;;
        -g|--hg) hg=$2; shift 2;;
        -b|--barcode) barcode=$2; shift 2;;
        -v|--vcf) vcf=$2; shift 2;;
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

sid=${smp_name}-${seq_type}-${phase_type}

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
raw_vname=${sid}.vcf.gz
raw_vpath=$vcf

if [ "$phase_type" == "phase" ]; then
    aim="filter SNPs: keep heterozygous SNPs only"
    flt_vname=${raw_vname/.vcf/.het.vcf}
    flt_vpath=$out_dir/$flt_vname
    cmd="$bin_bcftools view -Oz -i 'GT = \"het\"' $raw_vpath > $flt_vpath"
else
    aim="filter SNPs having low GP and keep heterozygous SNPs only"
    flt_vname=${raw_vname/.vcf/.gp.het.vcf}
    flt_vpath=$out_dir/$flt_vname
    cmd="$bin_bcftools view -Ou -i 'MAX(GP) > 0.99' $raw_vpath | 
      $bin_bcftools view -Oz -i 'GT = \"het\"' > $flt_vpath"
fi
eval_cmd "$cmd" "$aim"

aim="convert genome build, target build is hg$hg"
lift_vname=${flt_vname%.vcf.gz}.hg$hg.vcf.gz
lift_vpath=$out_dir/$lift_vname
if [ $hg -eq 19 ]; then
    log_msg "$aim"
    log_msg "Already hg19, skip liftover"
    lift_vname=$flt_vname
    lift_vpath=$flt_vpath
else
    cmd="$bin_python $bin_py_liftover -c $chain_hg19to38 -i $flt_vpath \\
           -o ${lift_vpath/.vcf/.tmp.vcf} -P $bin_liftover &&                    
         $bin_bcftools view -i 'POS > 0' -Oz ${lift_vpath/.vcf/.tmp.vcf} > ${lift_vpath} &&     
         rm ${lift_vpath/.vcf/.tmp.vcf}"
    eval_cmd "$cmd" "$aim"
fi

# chroms in Sanger Imputation Server genome reference have no leading 'chr'
aim="add leading chr for chrom names"
chr_vname=${lift_vname/.vcf/.chr.vcf}
chr_vpath=$out_dir/$chr_vname
cat $fasta | awk 'NR == 1 {print; exit}' | grep -i '^>chr'
if [ $? -eq 0 ]; then   # chrom names have leading chr
    cmd="$bin_bcftools annotate -Oz --rename-chrs $ensembl2ucsc $lift_vpath > $chr_vpath"
    eval_cmd "$cmd" "$aim"
else
    log_msg "$aim"
    log_msg "target fasta has no leading chr; skip this step"
    chr_vname=$lift_vname
    chr_vpath=$lift_vpath
fi

aim="bcftools fixref checking"
cmd="$bin_bcftools +fixref $chr_vpath -- -f $fasta"
eval_cmd "$cmd" "$aim"

aim="xcltk fixref"
fix_vname=${chr_vname/.vcf/.fixref.vcf}
fix_vpath=$out_dir/$fix_vname
tmp_prefix=${chr_vpath%.vcf.gz}
cmd="$bin_bcftools query -f '%CHROM:%POS-%POS\n' $chr_vpath > ${tmp_prefix}.region.lst &&  
     $bin_samtools faidx -r ${tmp_prefix}.region.lst $fasta | 
     $bin_bgzip -c > ${tmp_prefix}.fa.gz &&                     
     $bin_xcltk fixref -i $chr_vpath -r ${tmp_prefix}.fa.gz |  
     $bin_bgzip -c > $fix_vpath && 
     rm ${tmp_prefix}.region.lst"
eval_cmd "$cmd" "$aim"

aim="bcftools fixref checking"
cmd="$bin_bcftools +fixref $fix_vpath -- -f $fasta"
eval_cmd "$cmd" "$aim"

aim="filter duplicates (chrom + pos) and sort"
uniq_vname=${fix_vname/.vcf/.uniq.sort.vcf}
uniq_vpath=$out_dir/$uniq_vname
cmd="zcat $fix_vpath | awk '\$0 ~ /^#/ {print; next;} ! a[\$1\":\"\$2] {print; a[\$1\":\"\$2]=1}' |
     $bin_bcftools sort -Oz > $uniq_vpath"
eval_cmd "$cmd" "$aim"

csp_in_vname=$uniq_vname
csp_in_vpath=$uniq_vpath

aim="cellsnp pileup"
csp_dir=$out_dir/cellsnp
if [ "$seq_type" == "dna" ]; then
    cmd="$bin_cellsnp -s $bam -b $barcode -O $csp_dir -R $csp_in_vpath --minCOUNT 1 --minMAF 0 \\
      --minLEN 30 --minMAPQ 20 --inclFLAG 0 --exclFLAG 1796 --UMItag None -p $ncores  \\
      --genotype --gzip"
elif [ "$seq_type" == "rna" ]; then
    cmd="$bin_cellsnp -s $bam -b $barcode -O $csp_dir -R $csp_in_vpath --minCOUNT 1 --minMAF 0 \\
      --minLEN 30 --minMAPQ 20 --inclFLAG 0 --exclFLAG 772 --UMItag UR -p $ncores      \\
      --genotype --gzip"
else  # unknown
    log_msg "Warning: unknown seq type $seq_type, use the dna pileup method"
    cmd="$bin_cellsnp -s $bam -b $barcode -O $csp_dir -R $csp_in_vpath --minCOUNT 1 --minMAF 0  \\
      --minLEN 30 --minMAPQ 20 --inclFLAG 0 --exclFLAG 1796 --UMItag None -p $ncores           \\
      --genotype --gzip"
fi
eval_cmd "$cmd" "$aim"

csp_vpath=$csp_dir/cellSNP.base.vcf.gz

aim="merge pileup vcf and phase GT vcf"
gt_vname=${csp_in_vname/.vcf/.gt.vcf}
gt_vpath=$out_dir/$gt_vname
cmd="$bin_bcftools view -Oz -T $csp_vpath $csp_in_vpath > $gt_vpath"
eval_cmd "$cmd" "$aim"

#aim="extract phased GT"
#gt_tsv=$out_dir/${csp_in_vname%.vcf.gz}.gt.tsv.gz
#cmd="$bin_bcftools query -f '%CHROM\t%POS[\t%GT]\n' $gt_vpath | 
#     awk '{split(\$3, a, \"|\"); printf(\"%s\t%s\t%s\t%s\n\", \$1, \$2, a[1], a[2]); }' | 
#     $bin_bgzip -c > $gt_tsv"
#eval_cmd "$cmd" "$aim"

###### END ######
log_msg "All Done!"
log_msg "End"

