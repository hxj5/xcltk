#!/bin/awk -f
#get GQ score from GL.

{
    split($NF, a, ",");     # GL
    n = asort(a);
    if (n != 3) exit 1;

    gq = 10 * (a[3] - a[2])
    printf("%s\t%f\n", $0, gq);
}
