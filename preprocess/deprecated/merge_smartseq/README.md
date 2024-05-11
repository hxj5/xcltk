# Merge SMART-seq BAM Files

In [xcltk](https://github.com/hxj5/xcltk) preprocessing pipeline, the BAM 
files should be merged first if the input is SMART-seq data.

Here we give an example of how to merge the BAM files. It depends on some
softwares and files.

**Softwares**

- [samtools](https://github.com/samtools/samtools). 

**Files**

- A plain file listing pathes to all BAM files, one BAM file
per line. See 
[BCH869.input.492.bam.tsv](https://github.com/hxj5/xcltk/blob/master/preprocess/merge_smartseq/BCH869.input.492.bam.tsv) 
as an example.

Then merge the BAM files using `samtools merge`.
Note that `RG` tag should be added to make downstream analysis easier.
See details in the [samtools manual](http://www.htslib.org/doc/samtools-merge.html). 

```
# example: samtools merge -b BCH869.input.492.bam.tsv  -r  --write-index  -@ 3  BCH869.output.bam
samtools merge -b <input BAM list>  -r  --write-index  -@ 3  <output BAM>
```

One last thing is that the `barcode` file should be manually prepared, one barcode per cell. 
Specifically, the barcode should follow the string format specified by the `RG` tag 
in the output BAM file. 

For example, in the file `BCH869.output.bam`, each read is assigned a `RG` tag:

```
NS500400:224:H5YLJBGXY:3:13410:14263:15958      355     1       24591   1       38M     =       24723   170     CACTTCAGCCCCTCCCACACAGGGAAGCCAGATGGGTT  AAAAAEAEEEEEEEEE/EEEEEEEEEEEEEEEEEEEEE  NH:i:4  HI:i:2  AS:i:74 nM:i:0  RG:Z:BT_869-P07-A10.sort
```

Given the `RG` tag above (RG:Z:BT_869-P07-A10.sort), we see that a suffix `.sort` was added to the input file ID.
Following this rule, we can generate the `barcode` file based on the input BAM file list.
See [BCH869.output.492.RG.barcodes.tsv](https://github.com/hxj5/xcltk/blob/master/preprocess/merge_smartseq/BCH869.output.492.RG.barcodes.tsv) as an example.

Finally, the merged BAM file and manually generated barcode file would be used for downstream analysis.

