=======
History
=======

Release v0.1.11 (28/01/2021)
============================
| re-implement fixref with pysam:  
* support genome fasta as ref (-r)
* support gzip/bgzip input and output vcf
* support multiple alt alleles
* support multiple samples
* indels would be filtered
* support only ploidy = 2 for now

Release v0.1.10 (09/01/2021)
============================
* baf_post: support multiple BAMs
* baf_pileup: set cellTAG None when given bam list
* copy barcode file for baf_pileup and copy barcode & region
  files for phase_snp
* basefc: replace region.stop with region.end
* small fixes

Release v0.1.9 (04/01/2021)
===========================
* baf_pileup: add --uniqCOUNT
* specify sample ID through cmdline option

Release v0.1.8 (31/12/2020)
===========================
* phase_snp: fix load_phase
* baf_post: update pileup cmdline

Release v0.1.7 (29/12/2020)
===========================
* add pileup module and fix double counting

Release v0.1.6 (28/12/2020)
===========================
* phase_snp: support bed,gff,tsv for input region
* phase_snp: support vcf as input for phase file
* add gzip support for region sub-module
* baf_pre_impute: add -C/--call option and use cellsnp-lite
  by default to call germline SNPs instead of freebayes

Release v0.1.5 (19/12/2020)
===========================
* small fix
* baf_pre_impute and baf_pileup pass tests

Release v0.1.4 (17/12/2020)
===========================
* add baf_pileup pipeline

Release v0.1.3 (16/12/2020)
===========================
* add baf_pre_imputation pipeline

Release v0.1.2 (15/12/2020)
===========================
* add utils

Release v0.1.1 (14/12/2020)
===========================
* add fixref

Release v0.1.0 (13/12/2020)
===========================
* add feature-count

Release v0.0.2 (13/12/2020)
===========================
* add xcltk cmdline

Release v0.0.1 (12/12/2020)
===========================
* init modules: baf, rdr and reg
* add cmdline apps: xcltk-baf, xcltk-rdr and xcltk-reg
