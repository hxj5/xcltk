# TODO List

## Input


## Output


## BAF
### Add --minINCLUDE for read filtering


## RDR


## Implementation
### Use logging to replace `stdout.write()` and `stderr.write()`.


## Docs
Current xcltk feature counting (basefc) is close to htseq-count `union` mode 
with `--nonunique all`.
While counting each feature in parallel, xcltk basefc counts one read for both
features if the read effectively aligned to both features (i.e., passing all
read filtering).
Ref: https://htseq.readthedocs.io/en/release_0.11.1/count.html


## Scripts


## Discussion
### How to process some special UMIs in feature counting?
1. Within one UMI, reads overlap different features, e.g., due to errors in
   UMI collapsing, or features overlapping each other.
2. Within one UMI, some reads overlap one feature while others do not overlap
   any features.
3. Within one UMI, all reads overlap one feature, while some pass the --minINCLUDE 
   filtering and others do not.

