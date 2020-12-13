#!/bin/bash
#install all softwares needed for xcltk pipeline

if [ $# -lt 1 ]; then
    echo "Usage: $0 <all|pre|post>" >&2
    echo "Notes:" >&2
    echo "  all:  install all softwares of the pipeline" >&2
    echo "  pre:  install softwares needed only by the pre-imputation pipeline" >&2
    echo "  post: install softwares needed only by the post-imputation pipeline" >&2
    exit 1
fi

