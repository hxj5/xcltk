# TODO List

## Input


## Output


## BAF
### Add `xf` tag for filtering 10x reads.
### Add --minINCLUDE for read filtering


## RDR


## Implementation
### Add `xf` tag for filtering 10x reads.
### Use logging to replace `stdout.write()` and `stderr.write()`.


## Docs


## Scripts


## Discussion
### How to process some special UMIs in feature counting?
1. Within one UMI, reads overlap different features, e.g., due to errors in
   UMI collapsing, or features overlapping each other.
2. Within one UMI, some reads overlap one feature while others do not overlap
   any features.
3. Within one UMI, all reads overlap one feature, while some pass the --minINCLUDE 
   filtering and others do not.

