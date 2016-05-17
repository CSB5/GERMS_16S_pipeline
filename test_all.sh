#!/bin/bash

# http://redsymbol.net/articles/unofficial-bash-strict-mode/
set -euo pipefail

echo "WARN: make sure this is mostly python2" 1>&2

MYNAME=$(basename $(readlink -f $0))

usage() {
    echo "$MYNAME: run tests"
    echo " -d: Run dry-run tests"
    echo " -r: Run real-run tests"
}

skip_dry_runs=1
skip_real_runs=1
while getopts "dr" opt; do
    case $opt in
        d)
            skip_dry_runs=0
            ;;
        r)
            skip_real_runs=0
            ;;
        \?)
            usage
            exit 1
            ;;
    esac
done

args=""
if [ $skip_dry_runs -ne 1 ]; then
    args="$args -d"
fi
if [ $skip_real_runs -ne 1 ]; then
    args="$args -r"
fi
#echo "DEBUG args=$args" 1>&2

cd $(dirname $0)
commit=$(git describe --always --dirty)

echo "------------------------------------------------------------"
echo "Running static code checks with pylint in lib"
echo "------------------------------------------------------------"
set +e
for f in $(ls src/*py); do
    echo "Checking $f"
    PYTHONPATH=$(dirname $MYNAME)/lib pylint -j 2 -E --rcfile pylintrc $f
done
echo "Done"
set -e   


# e.g. odir=$(mktemp -d ../out/S16_XXXXXX); rmdir $odir; ./src/S16.py -1 ./data/p33_104k/p33_fwd.fastq.gz  -2 ./data/p33_104k/p33_rvs.fastq.gz -o $odir  -i 50 -m 100
echo "FIXME implement" 1>&2
