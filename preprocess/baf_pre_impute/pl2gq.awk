#!/bin/awk -f
#get GQ score from PL.

{
    split($NF, a, ",");     # PL
    n = asort(a);
    if (n != 3) exit 1;

    gq = a[2] - a[1]
    printf("%s\t%f\n", $0, gq);
}
