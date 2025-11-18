=======
History
=======


Release v0.5.2 (18/11/2025)
===========================
* baf: by default use ``min_count=11`` and ``min_maf=0.1`` for SNP filtering.



Release v0.5.1 (03/11/2025)
===========================
Fix bug (related to issue #13).

* baf.afc: filter SNPs in the input adata to match the input phased SNPs, 
  since the phasing step may filter some SNPs.
* utils.csp_io: transpose adata to save into standard feature x cell matrices.



Release v0.5.0 (17/09/2025)
===========================
This version mainly adds local phasing into the BAF module to fix the
errors of reference phasing in long genes.

Specifically, it modifies the local phasing function in XClone by updating
the probabilities of haplotype assignment in the EM iterations with 
Gaussian smoothing (default) or HMM smoothing along all SNPs within the gene,
considering the strong correlation between the haolotype assignment of closly 
adjacent SNPs (kind of to balance of prior biological knowledge and
data-driven evidence).

Notably, to balance reference phasing and local phasing, only a few long genes
(filtering with ``rlp_min_len``, ``rlp_min_n_snps``, and ``rlp_min_gap``)
will be considered for local phasing.

Additionally, preliminary analysis has shown that removing reference cells can
improve the local phasing results.
Generally, gene BAFs in reference cells are near 0.5, which has little 
information about local phasing and can be treated as noise.
A new optional parameter ``ref_cell_fn`` is added to specify the barcodes
of reference cells.


BAF module

* post-genotype filtering SNPs given ``min_count`` and ``min_maf``: 
  considering only REF and ALT (AD & DP) alleles, but not OTH alleles.
  Cellsnp-lite (at least v1.2.3 and before) may keep some SNPs unexpectly,
  e.g., SNP with ``AD=9;DP=9;OTH=1`` when ``minMAF=0.1; minCOUNT=10``.
* baf.pipeline: add options ``--minCOUNT`` and ``--minMAF``.
* afc: sort SNP list in ``reg.snp_list``.
* baf: delete afc cmdline interface and remove ``afc`` and ``rpc`` from 
  xcltk subcommand.
* restructure: move the ``baf.count`` to ``baf.fc.main``.


RDR module

* restructure: move the ``rdr.count`` to ``rdr.fc.main``.


Others

* utils: add sub-module ``csp_io`` for loading and saving cellsnp-lite outputs.
* restructure: move the ``rdr.count`` to ``rdr.fc.main``.
* preprocess: update file download links.
* docs: add brief summary of xcltk basefc and link it to the ``htseq-count``.
* docs: update TODO list.
* README: add citation.



Release v0.4.1 (10/01/2025)
===========================
* RDR: `--minINCLUDE` supports both INT and FLOAT.
  If INT, the minimum length of included part within specific feature.
  If FLOAT, the minimum fraction of included part within specific feature.
* BAF: add index to the output folder of each step, e.g., "1_pileup".
* docs: update TODO list. 


Release v0.4.0 (06/12/2024)
===========================
* RDR: add `--minINCLUDE` option for read filtering, which is the minimum 
  length of included part within specific feature.
  For example, if the genomic range of a feature is chr1:1000-3000, and one
  fetched read (100bp) aligned to two locus, chr1:601-660 (60bp) and
  chr1:3801-3840 (40bp), then no any part of the read is actually included
  within the feature, hence it will be filtered by `--minINCLUDE=30`, whereas
  older versions of xcltk may keep the read.
  Note, as the feature counting in RDR is performmed independently for each
  feature, so one read filtered by `--minINCLUDE` in one feature may still be
  fetched and counted by other features.
* update docstring, using the numpydoc style.
* add TODO list in `docs/TODO.md`.
* fix typo.


Release v0.3.1 (06/06/2024)
===========================
* BAF: in ref_phasing, use multiprocessing to phase SNPs of one chromosome
  per subprocess.
* specify dtype of column 0 as str in pd.read_csv() when loading region file.


Release v0.3.0 (11/05/2024)
===========================
The v0.2.x was skipped since this new version has several substantial updates:

* BAF: do reference phasing on local machines instead of using online 
  service.
* BAF & RDR: better support well-based (e.g., SMART-seq) data without
  the need to merge the input BAM files first;
* coding improvement using a more unified framework, mainly using the
  ``fc`` (feature counting) and ``utils`` sub-modules.

Feature enhancement

BAF part:

* add ``xcltk baf`` command line tool to support reference phasing on
  local machines instead of using online service.
* ``xcltk allelefc``: better support well-based (e.g., SMART-seq) data without
  the need to merge the input BAM files first;
* ``xcltk allelefc``: both REF and ALT allele counting will exclude the 
  UMIs/reads mapped to both alleles when *no_dup_hap* is True.

RDR part:

* better support well-based (e.g., SMART-seq) data without the need to merge
  the input BAM files first;
* re-implement the ``xcltk basefc`` using the ``fc`` (feature counting) 
  framework.

Preprocess:

* re-implement the preprocess pipeline by 
  (1) replace the bash scripts with python functions, e.g., 
      wrapping SNP calling (previously ``baf_pre_phase.sh``) into 
      ``xcltk.baf.genotype::pileup()``; 
      reference phasing locally with ``xcltk.baf.genotype::ref_phasing()``;
      wrapping allele-specific feature counting (previously 
      ``baf_post_phase.sh``) with ``xcltk.baf.count::afc_wrapper()``.
  (2) further wrap the three functions into a pipeline implemented as
      a sub-module ``xcltk.baf.pipeline`` and also as a command line tool
      ``xcltk baf``.

Others:

* rename the cmdline command ``xcltk pileup`` to ``xcltk allelefc``.
* make the cmdline options more unified, e.g., "--samList" and "--ncores" in
  "xcltk allelefc", "xcltk pipeline", and "xcltk basefc".
* usage() functions by default output to *stdout* instead of *stderr*.
* cmdline "--help" option exit code changes from 1 to 0.
* add/update a few util sub-modules such as ``vcf.py``, ``xlog.py`` etc.
* add post_hoc scripts for post-processing xcltk output.
* initialize "data" dir and add feature annotation files.


Release v0.1.16 (28/01/2023)
============================
* baf: add reference phasing correction (xcltk rpc).
* preprocess: restructure, update scripts and data.
* rdr: output 4-column features.

Release v0.1.15 (17/07/2022)
============================
* baf_pre_impute: keep het SNPs only after calling germline SNPs
* baf_post_impute: output all regions when running xcltk pileup
* rdr: fix a bug that pysam was not imported.

Release v0.1.14 (01/05/2022)
============================
| update baf haploblock pileup:
* re-implement the module
* fix the double counting issue of UMIs or reads when aggregating phased SNPs
  (some UMIs or reads could cover more than one SNPs)
* fix the issue that some UMIs are aligned to both haplotype alleles
  (--countDupHap)
* add an option to output all regions (--outputAllReg)

Release v0.1.13 (02/08/2021)
============================
* rdr: fix program suspension caused by unmatched chrom

Release v0.1.12 (26/02/2021)
============================
* baf_pre: add --umi and --duplicates options

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
