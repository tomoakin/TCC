#!/bin/sh
export TCC_REAL_TEST=true
for f in TCC/tests/*.R
do
    echo ""
    echo "==== Running $f"
    time R --slave < $f
done

echo "Done"
