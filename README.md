# xcltk: Toolkit for XClone Preprocessing

[![](https://img.shields.io/pypi/v/xcltk.svg)][pypi]
[![](https://img.shields.io/github/license/hxj5/xcltk)][licence]


[XClone][XClone repo] is a statistical method to detect allele- and 
haplotype-specific copy number variations (CNVs) and reconstruct 
tumour clonal substructure from scRNA-seq data, 
by integrating the expression levels (read depth ratio; RDR signals) and 
the allelic balance (B-allele frequency; BAF signals).
It takes three matrices as input: the allele-specific *AD* and *DP* matrices
(BAF signals) and the *total read depth* matrix (RDR signals).

The [xcltk][xcltk repo] package implements a preprocessing pipeline to 
generate the three matrices from SAM/BAM/CRAM files.
It supports data from multiple single-cell sequencing platforms, including 
droplet-based (e.g., 10x Genomics) and well-based (e.g., SMART-seq)
platforms.

For details of xcltk and XClone, please checkout our paper:

> Huang, R., Huang, X., Tong, Y. et al. Robust analysis of allele-specific copy number alterations from scRNA-seq data with XClone. Nat Commun 15, 6684 (2024). https://doi.org/10.1038/s41467-024-51026-0


## News

You can find the full manual of the xcltk preprocessing pipeline at
[preprocess/README.md][preprocess manual].

All release notes are available at [docs/release.rst][release]


## Installation

### Install via pip (latest stable version)

xcltk is avaliable through [pypi][pypi].

```shell
pip install -U xcltk
```

### Install from this Github Repo (latest stable/dev version)

```shell
pip install -U git+https://github.com/hxj5/xcltk
```

In either case, if you don't have write permission for your current Python
environment, we suggest creating a separate [conda][conda] environment 
or add `--user` for your current one.


## Manual

You can check the full parameters with `xcltk -h`.

```
Program: xcltk (Toolkit for XClone Preprocessing)
Version: 0.5.2

Usage:   xcltk <command> [options]

Commands:
  -- BAF calculation
     baf              Preprocessing pipeline for XClone BAF.
     fixref           Fix REF allele mismatches based on reference FASTA.

  -- RDR calculation
     basefc           Basic feature counting.

  -- Tools
     convert          Convert between different formats of genomic features.

  -- Others
     -h, --help       Print this message and exit.
     -V, --version    Print version and exit.
```



[conda]: https://docs.conda.io/en/latest/
[licence]: https://github.com/hxj5/xcltk
[preprocess manual]: https://github.com/hxj5/xcltk/tree/master/preprocess
[pypi]: https://pypi.org/project/xcltk
[release]: https://github.com/hxj5/xcltk/blob/master/docs/release.rst
[XClone repo]: https://github.com/single-cell-genetics/XClone
[xcltk repo]: https://github.com/hxj5/xcltk
