# utils.sh - Utils
# Author: Xianjie Huang


# check if file exists
# @example assert_e  "$file"  "Test File"
function assert_e() {
    if [ ! -e "$1" ]; then
        echo "$2 $1 does not exist!" >&2
        exit 1
    fi
}


# check if string is not empty
# @example assert_n  "$sid"  "Sample ID"
function assert_n() {
    if [ -z "$1" ]; then
        echo "$2 $1 is empty!" >&2
        exit 1
    fi
}


function generate_bin_gl2gq() {
    cat <<EOF > $1
#!/bin/awk -f
#convert GL to GQ score.
{
    split(\$NF, a, ",");     # GL
    n = asort(a);
    if (n != 3) exit 1;

    gq = 10 * (a[3] - a[2])
    printf("%s\t%f\n", \$0, gq);
}
EOF
    chmod u+x $1
}


function generate_bin_pl2gq() {
    cat <<EOF > $1
#!/bin/awk -f
#convert PL to GQ score.
{
    split(\$NF, a, ",");     # PL
    n = asort(a);
    if (n != 3) exit 1;

    gq = a[2] - a[1]
    printf("%s\t%f\n", \$0, gq);
}
EOF
    chmod u+x $1
}


function __now() {
    date '+%Y-%m-%d %H:%M:%S'
}


function __exit() {
    echo "[`__now`] $1" >&2
    exit 1
}


function log_err() {
    echo "[`__now`] $1" >&2
}


function log_msg() {
    echo "[`__now`] $1"
}

