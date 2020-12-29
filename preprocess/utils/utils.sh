# Utils
# Author: Xianjie Huang

function __now() {
    date '+%Y-%m-%d %H:%M:%S'
}

function __exit() {
    echo "[`__now`] $1" >&2
    exit 1
}

#@abstract  Execute command and check running status
#@param $1  Command [str]
#@param $2  Aim of this command [str]
#@return    No RetCode
#@note      Will exit the program if error.
#@example   eval_cmd "echo hello world" "test this function"
function eval_cmd() {
    if [ $# -lt 2 ]; then
        __exit "[E::eval_cmd] too few arguments."
    fi
    local cmd="$1"
    local aim="$2"

    log_msg "START $aim"
    echo "=> COMMAND"
    echo "$cmd"
    echo "=> OUTPUT"
    eval "$cmd"
    echo "=> DONE"
    if [ $? -ne 0 ]; then
        __exit "[E::eval_cmd] failed to run PART $aim"
    fi
    log_msg "END $aim"
    echo
}

function load_cfg() {
    if [ $# -lt 1 ]; then
        __exit "[E::load_cfg] too few arguments."
    fi
    local cfg=$1
    local cmd=
    while read line; do
        cmd=`echo "$line" | sed 's/ *= */=/'`
        eval "$cmd"
    done < $cfg
}

function log_err() {
    echo "[`__now`] $1" >&2
}

function log_msg() {
    echo "[`__now`] $1"
}

