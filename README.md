# xcltk: Toolkit for XClone Preprocessing

[![](https://img.shields.io/pypi/v/xcltk.svg)][pypi]
[![](https://img.shields.io/github/license/hxj5/xcltk)][licence]


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


## News

You can find the full manual of the xcltk preprocessing pipeline at
[preprocess/README.md][preprocess manual].

All release notes are available at [docs/release.rst][release]


## Installation

xcltk is avaliable through [pypi][pypi].
To install, type the following command line, and add `-U` for upgrading:

```shell
pip install -U xcltk
```

Alternatively, you can install from this GitHub repository for latest (often
development) version by following command line

```shell
pip install -U git+https://github.com/hxj5/xcltk
```

In either case, if you don't have write permission for your current Python 
environment, we suggest creating a separate [conda][conda] environment 
or add --user for your current one.


## Manual

You can check the full parameters with `xcltk -h`.

```
Program: xcltk (Toolkit for XClone Preprocessing)
Version: 0.3.0

Usage:   xcltk <command> [options]

Commands:
  -- BAF calculation
     allelefc         Allele counting for each feature.
     baf              Preprocessing pipeline for XClone BAF.
     fixref           Fix REF, ALT and GT.
     rpc              Reference phasing correction.

  -- RDR calculation
     basefc           Basic feature counting.

  -- Tools
     convert          Convert different region file formats.

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