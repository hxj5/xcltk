# xcltk: Toolkit for XClone Preprocessing

[XClone][XClone repo] is a statistical method to detect haplotype-aware 
copy number variations (CNVs) and reconstruct tumour clonal substructure from
scRNA-seq data, 
by integrating the expression levels (read depth ratio; RDR signals) and 
the allelic balance (B-allele frequency; BAF signals).
It takes three matrices as input: the allele-specific *AD* and *DP* matrices
(BAF signals) and the *total read depth* matrix (RDR signals).

The [xcltk][xcltk repo] package implements a preprocessing pipeline to 
generate the three matrices from SAM/BAM/CRAM files.
It supports data from multiple single-cell sequencing platforms, including 
droplet-based (e.g., 10x Genomics) and well-based (e.g., SMART-seq)
platforms.


## Installation

### Softwares

To use the pipeline, please install
[python (tested on python 3.11)](https://www.python.org/) and 
[xcltk >= 0.3.0][xcltk repo], together with a few dependencies listed below.

- [bcftools][bcftools]
- [bgzip or htslib][htslib]
- [cellsnp-lite >= 1.2.0][cellsnp-lite]
- [eagle2][eagle2]

Please install these softwares and add them to the system search path (i.e.,
the system variable PATH).

#### Use conda environment (recommended)

To make the installation easier, we highly recommend to set up a conda env and
then install the softwares in the env.

```
# Notes:
# 1. `pandas` supports python 3.9-3.12 while `pysam` supports python 3.6-3.11,
#     therefore, we tested with python 3.11;
# 2. bgzip/htslib will be automatically installed when installing bcftools
#    via conda.
# 3. Eagle2 has to be manually installed since there is not corresponding
#    installation in conda.

conda create -n xcltk python=3.11
conda activate xcltk
conda install -c conda-forge -c bioconda bcftools cellsnp-lite
pip install 'xcltk>=0.3.0'
```

Importantly, [eagle2][eagle2] has to be manually installed since there is not 
corresponding package in conda.
You can download the
[Eagle software](https://alkesgroup.broadinstitute.org/Eagle/#x1-30002)
and then unzip the file, e.g.,

```
# Commands below download and unzip Eagle v2.4.1;
# Download the latest version at https://alkesgroup.broadinstitute.org/Eagle/

wget https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/Eagle_v2.4.1.tar.gz
tar xzvf Eagle_v2.4.1.tar.gz
```

### Files

In addition, the pipeline relies on a common SNP VCF (for pileup), 
phasing reference panel and genetic map (for phasing with eagle2).
You can download the files using following links.

#### 1000G SNP VCF

```
# hg38
wget https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz

# hg19
wget https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e2.chr1toX.hg19.vcf.gz
```

#### 1000G Reference Panel

The pre-compiled files below will be used by Eagle2 as reference panel for
SNP phasing.
Credits to the authors of [Numbat][Numbat].

```
# hg38
wget http://pklab.med.harvard.edu/teng/data/1000G_hg38.zip

# hg19
wget http://pklab.med.harvard.edu/teng/data/1000G_hg19.zip
```

#### Genetic Map

Use commands in above section to download zipped Eagle2 file.
After unzip, the genetic map files are in subfolder `tables`,
e.g., `Eagle_v2.4.1/tables`.

#### Feature annotations

Check [data/anno](https://github.com/hxj5/xcltk/tree/master/data/anno) folder
in the [xcltk][xcltk repo] Github repo for the latest version of feature 
annotation files.

```
# hg38
wget https://raw.githubusercontent.com/hxj5/xcltk/master/data/anno/annotate_genes_hg38_update_20230126.txt

# hg19
wget https://raw.githubusercontent.com/hxj5/xcltk/master/data/anno/annotate_genes_hg19_update_20230126.txt
```


## Quick Start

As mentioned, the [xcltk][xcltk repo] preprocessing pipeline aims to generate
the three matrices, the *AD*, *DP*, and *total read depth* matrices, 
from input BAM file(s).
The allele-specific AD and DP matrices, refered to as **BAF** 
(B-Allele Frequency), contains the information of allelic balance.
The total read depth, refered to as **RDR** (read depth ratio), 
contains the information of gene expression levels.

For BAF part, xcltk implements a `xcltk baf` command line tool,
while for RDR part, you may use any available feature counting tools, 
e.g., [cellranger][cellranger] and [STARsolo][STARsolo], or simply
use `xcltk basefc` command line tool.


### 10x scRNA-seq data

#### BAF part

An example for **hg38** data.
Type `xcltk baf --help` for full parameters.

```
# conda activate xcltk

xcltk baf  \
    --label        {sample name}    \
    --sam          {BAM file}       \
    --barcode      {barcode file}   \
    --snpvcf       {genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz}  \
    --region       {annotate_genes_hg38_update_20230126.txt}        \
    --outdir       {output folder}          \
    --gmap         {Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz}  \
    --eagle        {Eagle_v2.4.1/eagle}      \
    --paneldir     {1000G_hg38}              \
    --ncores       10
```

#### RDR part

An example for **hg38** data.
Type `xcltk basefc --help` for full parameters.

```
# conda activate xcltk

xcltk basefc     \
    --sam          {BAM file}       \
    --barcode      {barcode file}   \
    --region       {annotate_genes_hg38_update_20230126.txt}   \
    --outdir       {output folder}   \
    --ncores       10
```


### SMARTA-seq data

#### BAF part

An example for **hg19** data.
Type `xcltk baf --help` for full parameters.

```
# conda activate xcltk

xcltk baf  \
    --label        {sample name}           \
    --samList      {BAM list file}         \
    --sampleList   {sample ID list file}   \
    --snpvcf       {genome1K.phase3.SNP_AF5e2.chr1toX.hg19.vcf.gz}  \
    --region       {annotate_genes_hg19_update_20230126.txt}        \
    --outdir       {output folder}          \
    --gmap         {Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz}  \
    --eagle        {Eagle_v2.4.1/eagle}      \
    --paneldir     {1000G_hg19}              \
    --cellTAG      None         \
    --UMItag       None         \
    --ncores       10
```

#### RDR part

An example for **hg19** data.
Type `xcltk basefc --help` for full parameters.

```
# conda activate xcltk

xcltk basefc     \
    --samList      {BAM list file}       \
    --sampleList   {sample ID list file}   \
    --region       {annotate_genes_hg19_update_20230126.txt}   \
    --outdir       {output folder}   \
    --cellTAG      None         \
    --UMItag       None         \    
    --ncores       10
```


### Notes

As elaborated in the next section, the BAF calculation is mainly implemented
in three functions:

- `pileup()` function in module `xcltk.baf.genotype` for germline heterozygous
  SNP calling.
- `ref_phasing()` function in module `xcltk.baf.genotype` for SNP
  reference phasing with Eagle2.
- `afc_wrapper()` function in module `xcltk.baf.count` for allele-specific
  feature counting.

The `xcltk baf` command line is actually a wrapper of these three functions.
You may select part of the three functions and update them on your needs
(refer to the 
[xcltk.baf.pipeline](https://github.com/hxj5/xcltk/blob/master/xcltk/baf/pipeline.py)
script).

**When matched omics data is available** 
If you have matched omics data, especially the DNA data, you may try using
it to get the phased SNPs.
For example, if you have matched scDNA-seq data of the same sample,
you may

1. call `xcltk.baf.genotype::pileup()` function on the scDNA-seq data 
   to obtain a more reliable genotyping results.
2. call `xcltk.baf.genotype::ref_phasing()` function to do the reference
   phasing.
   Before phasing, you may need to split the VCF file by chromosomes 
   with `xcltk.utils.vcf::vcf_split_chrom()` and index the new VCF files with
   `xcltk.utils.vcf::vcf_index()`.
   After phasing, you may merge the chromosome-specific phased VCF files with 
   `xcltk.utils.vcf::vcf_merge()`.
3. call `xcltk.baf.count::afc_wrapper()` function to do allele-specific
   feature counting on the scRNA-seq data.


## Pipeline walkthrough

### BAF part

[xcltk][xcltk repo] implements `xcltk baf` command line tool to compute 
the BAF signals for each feature and single cell.
It mainly includes three steps, while each step has been implemented in a
specific function.

#### 1. Germline SNPs calling

xcltk calls germline heterozygous SNPs using a pseudo-bulk strategy, i.e.,
by aggregating UMIs/reads of all cells as one bulk sample.
By default, it only keeps SNPs with `--minCOUNT 20` and `--minMAF 0.1`.

This step is implemented in the function `pileup()` in module
`xcltk.baf.genotype`.

#### 2. Reference based SNP phasing

xcltk performs SNP phasing with Eagle2 using 1000G reference panel.

This step is implemented in the function `ref_phasing()` in module
`xcltk.baf.genotype`.

#### 3. BAF calculation for features and single cells

xcltk computes the BAF signals by aggregating the counts of unique
haplotype-specific UMIs or reads, which are fetched from phased SNPs, 
in haplotype features (blocks) and single cells.

This step is implemented in the function `afc_wrapper()` in module
`xcltk.baf.count` and also as a command line tool `xcltk allelefc`.

Finally, matrices of *feature x cell* AD (UMI or read counts of one haplotype)
and DP (UMI or read counts of both haplotypes) would be outputed
for downstream analysis.

### RDR part

[xcltk][xcltk repo] computes the RDR signals in each feature and cell by 
extracting and counting UMIs or reads overlapping with target features, 
during which minimum reads filtering was performed.

By default, reads with MAPQ<20, aligned length<30nt are filtered. 
In addition, reads with certain FLAG bits will also be filtered.
Specifically, by default, reads with `UNMAP,SECONDARY,QCFAIL` 
(when using UMI; e.g., 10x scRNA-seq data) or 
`UNMAP,SECONDARY,QCFAIL,DUP` (when not using UMIs; e.g., SMART-seq2, 
10x scDNA-seq, 10x scATAC-seq) will be filtered.
Intronic reads and reads that are mapped to multiple features are largely 
preserved for counting. 
Finally, a *feature x cell* matrix of UMI or read counts would be outputed 
for downstream analysis.

This step is implemented in the `xcltk basefc` command line tool.



[bcftools]: https://github.com/samtools/bcftools
[cellranger]: https://www.10xgenomics.com/support/software/cell-ranger/latest
[cellsnp-lite]: https://github.com/single-cell-genetics/cellsnp-lite
[eagle2]: https://alkesgroup.broadinstitute.org/Eagle/
[htslib]: https://github.com/samtools/htslib
[Numbat]: https://github.com/kharchenkolab/numbat/
[STARsolo]: https://github.com/alexdobin/STAR
[XClone repo]: https://github.com/single-cell-genetics/XClone
[xcltk repo]: https://github.com/hxj5/xcltk