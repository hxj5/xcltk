# xcltk

[![](https://img.shields.io/pypi/v/xcltk.svg)][pypi]
[![](https://img.shields.io/github/license/hxj5/xcltk)][licence]

xcltk is a toolkit to preprocess the BAF and RDR signals from single cell 
data for CNV detection by [XClone][XClone].
You can find the scripts of preprocessing at [preprocess][preprocess].

All release notes are available at [docs/release.rst][release]

## Installation

xcltk is avaliable through [pypi][pypi]. To install, type the following command 
line, and add `-U` for upgrading:

```shell
pip install -U xcltk
```

Alternatively, you can install from this GitHub repository for latest (often 
development) version by following command line

```shell
pip install -U git+https://github.com/hxj5/xcltk
```

In either case, if you don't have write permission for your current Python environment,
 we suggest creating a separate [conda][conda] environment or add --user for your 
current one.

## Manual

You can check the full parameters with `xcltk -h`.

```
Program: xcltk (Toolkit for XClone)
Version: 0.1.16

Usage:   xcltk <command> [options]

Commands:
  -- BAF calculation
     fixref           Fix REF, ALT and GT.
     pileup           Pileup, support unique counting.
     rpc              Reference phasing correction.

  -- RDR calculation
     basefc           Basic feature counting.

  -- Region operations
     convert          Convert different region file formats.

  -- Others
     -h, --help       Print this message and exit.
     -V, --version    Print version and exit.
```

[pypi]: https://pypi.org/project/xcltk
[licence]: https://github.com/hxj5/xcltk
[release]: https://github.com/hxj5/xcltk/blob/master/doc/release.rst
[conda]: https://docs.conda.io/en/latest/
[XClone]: https://github.com/single-cell-genetics/XClone
[preprocess]: https://github.com/hxj5/xcltk/tree/master/preprocess

