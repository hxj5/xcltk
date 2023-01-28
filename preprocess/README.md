# XClone Preprocessing

[XClone][XClone repo] is a statistical method that 
integrates expression levels (RDR signals) and allelic balance (BAF signals) to 
detect haplotype-aware CNVs and 
reconstruct tumour clonal substructure from scRNA-seq data. It takes three matrices 
as input: allele-specific AD and DP matrices (BAF signals) and the total read depth 
matrix (RDR signals).

This preprocessing pipeline is aimed to generate the three 
matrices from SAM/BAM/CRAM files. It supports data from multiple single-cell 
sequencing platforms, including 10x scRNA-seq and SMART-seq.

## Dependencies

We have wrapped the preprocessing pipeline into a github repo
[xcltk][xcltk repo]
To use the pipeline, the repo should be cloned firstly.

```shell
git clone https://github.com/hxj5/xcltk.git
```

You can find the scripts and data for preprocessing in the directory 
[xcltk/preprocess][preprocess dir].

In addition, the pipeline dependends on the below softwares and files.

### Softwares

All tools below, after being installed, should be added to system search path (i.e.,
the system variable PATH).

- [bcftools](https://github.com/samtools/bcftools)
- [bgzip or htslib](https://github.com/samtools/htslib)
- [cellsnp-lite](https://github.com/single-cell-genetics/cellsnp-lite)
- [LiftOver](https://genome-store.ucsc.edu/)
- [python](https://www.python.org/)
- [samtools](https://github.com/samtools/samtools)
- [xcltk][xcltk repo]

Note that except `liftOver`, all tools above can be easily installed through `conda`
or `pip`.

### Files

[xcltk](https://github.com/hxj5/xcltk) uses 
[Sanger Imputation Server][Sanger Server] for 
reference phasing, i.e., infering haplotypes given a reference population panel.
For now, the server only supports hg19 and it provides a specific hg19 fasta file
listed below:

- [hg19 fasta](https://imputation.sanger.ac.uk/?resources=1) from Sanger Imputation 
  Server. Download from link 
  ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz

## Quick Start

The [xcltk][xcltk repo] preprocessing pipeline includes 2 parts: the RDR part and 
the BAF part.
The total read depth, refered to as RDR, contains the information of gene 
expression levels. 
The allele-specific AD and DP matrices, refered to as BAF (B-Allele Frequency), 
contains the information of allelic balance.

### 10x scRNA-seq data

#### RDR part

We recommend using `xcltk basefc` to generate the RDR matrix. It takes the 10x 
scRNA-seq data as input and will output a cell x gene UMI count sparse matrix.

Below is an example for 10x scRNA-seq data:

```
xcltk basefc               \
    -s  <BAM file>          \
    -b  <barcode file>      \
    -r  <region file>       \
    -T  tsv                 \
    -O  <out dir>           \
    -p  10                  \
    --cellTAG  CB           \
    --UMItag   UB           \
    --minLEN   30           \
    --minMAPQ  20           \
    --maxFLAG  4096
```

Alternatively, you can use the cell x gene UMI count matrix outputed by 
any other counting method, such as CellRanger.

#### BAF part

The BAF part includes 3 steps: pre-phasing (pre-Sanger-Phasing), Sanger Phasing,
and post-phasing (post-Sanger-Phasing). As mentioned above, xcltk uses
[Sanger Phasing Server][Sanger Server] for reference phasing.

##### Pre-Phasing

This step is aimed to call germline SNPs from input BAM file and save the SNPs
in a VCF file. To run this step, simply use the script
[baf_pre_phase.sh][pre-phasing script]. 

Below is an example for 10x scRNA-seq data, assuming it is aligned to hg38:

```
./baf_pre_phase.sh            \
    -N  <sample name>         \
    -s  <BAM file>            \
    -b  <barcode file>        \
    -F  <sanger hg19 fasta>   \
    -O  <out dir>             \
    -g  38                    \
    -C  CB                    \
    -u  UB                    \
    -p  10 
```

Run `./baf_pre_phase.sh -h` to check detailed information for each option.

#### Sanger Phasing

The VCF generated from previous step would be used for reference phasing on
[Sanger Imputation Server][Sanger Server]. There's a comprehensive introduction to
this step at [here][Sanger Wiki].

Briefly, The heterozygous SNPs, submitted to the server in a VCF file, are phased 
by selecting the "Haplotype Reference Consortium (r1.1)" reference panel and the 
"phase with EAGLE2, no imputation" pipeline.

The server will return a VCF file containing all phased SNPs.

#### Post-Phasing

This step is aimed to pileup UMI counts (AD and DP) for each SNP and then aggregate
UMI counts in gene level, using the Sanger Phasing results, to obtain allele-specific
cell x gene AD and DP matrices.

To run this step, simply use the script [baf_post_phase.sh][post-phasing script]

Below is an example for 10x scRNA-seq data, assuming it is aligned to hg38:

```
./baf_post_phase.sh        \
    -N  <sample name>      \
    -s  <BAM file>         \
    -b  <barcode file>     \
    -v  <phased VCF>       \
    -f  <fasta file>       \
    -O  <out dir>          \
    -g  38                 \
    -C  CB                 \
    -u  UB                 \
    -p  10
```

Run `./baf_post_phase.sh -h` to check detailed information for each option.

### SMARTA-seq data

As [xcltk][xcltk repo] only supports single BAM file as input for now,
the BAM files should be merged into one single BAM file first if the input is
SMART-seq data. Please follow the 
[example](https://github.com/hxj5/xcltk/tree/master/preprocess/merge_smartseq)
to merge the SMART-seq BAM files with proper settings.

Afterwards, there should be a merged BAM file containing `RG` tag and 
a manually generated barcode file.

#### RDR part

We recommend using `xcltk basefc` to generate the RDR matrix. It takes the merged 
BAM as input and will output a cell x gene read count matrix.

Below is an example for merged SMART-seq data:

```
xcltk basefc               \
    -s  <BAM file>          \
    -b  <barcode file>      \
    -r  <region file>       \
    -T  tsv                 \
    -O  <out dir>           \
    -p  10                  \
    --cellTAG  RG           \
    --UMItag   None         \
    --minLEN   30           \
    --minMAPQ  20           \
    --maxFLAG  255
```

Alternatively, you can use the cell x gene UMI count matrix outputed by 
any other counting method.

#### BAF part

The BAF part includes 3 steps: pre-phasing (pre-Sanger-Phasing), Sanger Phasing,
and post-phasing (post-Sanger-Phasing). As mentioned above, xcltk uses
[Sanger Phasing Server][Sanger Server] for reference phasing.

##### Pre-Phasing

This step is aimed to call germline SNPs from input BAM file and save the SNPs
in a VCF file. To run this step, simply use the script
[baf_pre_phase.sh][pre-phasing script]. 

Below is an example for merged SMART-seq data, assuming it is aligned to hg19:

```
./baf_pre_phase.sh            \
    -N  <sample name>         \
    -s  <BAM file>            \
    -b  <barcode file>        \
    -F  <sanger hg19 fasta>   \
    -O  <out dir>             \
    -g  19                    \
    -C  RG                    \
    -u  None                  \
    -p  10                    \
    --noDUP                   \
    --smartseq
```

Run `./baf_pre_phase.sh -h` to check detailed information for each option.

#### Sanger Phasing

The VCF generated from previous step would be used for reference phasing on
[Sanger Imputation Server][Sanger Server]. There's a comprehensive introduction to
this step at [here][Sanger Wiki].

Briefly, The heterozygous SNPs, submitted to the server in a VCF file, are phased 
by selecting the "Haplotype Reference Consortium (r1.1)" reference panel and the 
"phase with EAGLE2, no imputation" pipeline.

The server will return a VCF file containing all phased SNPs.

#### Post-Phasing

This step is aimed to pileup read counts (AD and DP) for each SNP and then aggregate
read counts in gene level, using the Sanger Phasing results, to obtain allele-specific
cell x gene AD and DP matrices.

To run this step, simply use the script [baf_post_phase.sh][post-phasing script]

Below is an example for 10x scRNA-seq data, assuming it is aligned to hg19:

```
./baf_post_phase.sh        \
    -N  <sample name>      \
    -s  <BAM file>         \
    -b  <barcode file>     \
    -v  <phased VCF>       \
    -f  <fasta file>       \
    -O  <out dir>          \
    -g  19                 \
    -C  RG                 \
    -u  None               \
    -p  10                 \
    --noDUP
```

Run `./baf_post_phase.sh -h` to check detailed information for each option.

### Notes

The scripts in the BAF part and RDR part can be modified to be applied on other
omics data, such as scDNA-seq or scATAC-seq data.

For example, if you have matched scDNA-seq data with the scRNA-seq data, it is
recommended to run the `Pre-Phasing` step using scDNA-seq data to obtain a
more reliable genotyping results.

## Details

### BAF part

[xcltk][xcltk repo] computes the BAF in each haplotype block and cell in three steps: 

#### 1. Germline SNPs calling

xcltk calls heterozygous germline SNPs by treating all cells as one bulk sample. 
UMI information (if avaiable) would be ignored in this step.

#### 2. Reference based SNP phasing

xcltk performs reference based SNP phasing by using Sanger Imputation Service 
(https://imputation.sanger.ac.uk/). The heterozygous SNPs,
submitted to the server in a VCF file, are phased by selecting the 
"Haplotype Reference Consortium (r1.1)" reference panel and the 
"phase with EAGLE2, no imputation" pipeline.

#### 3. BAF calculation in haplotype blocks and single cells

xcltk computes the BAF by aggregating the counts of **unique** haplotype-specific 
UMIs or reads, which are fetched from phased SNPs, 
in haplotype blocks and single cells. 

Finally, matrices of "block x cell" AD (UMI or read counts of one haplotype) and 
DP (UMI or read counts of both haplotypes) would be outputed
for downstream analysis.

These three steps are implemented and wrapped into several script files, 
which are publicly available at [xcltk/preprocess dir][preprocess dir].


[post-phasing script]: https://github.com/hxj5/xcltk/blob/master/preprocess/baf_post_phase.sh
[pre-phasing script]: https://github.com/hxj5/xcltk/blob/master/preprocess/baf_pre_phase.sh
[preprocess dir]: https://github.com/hxj5/xcltk/tree/master/preprocess
[Sanger Server]: https://imputation.sanger.ac.uk/
[Sanger Wiki]: https://imputation.sanger.ac.uk/?instructions=1
[XClone repo]: https://github.com/single-cell-genetics/XClone
[xcltk repo]: https://github.com/hxj5/xcltk

