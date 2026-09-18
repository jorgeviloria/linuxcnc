#!/bin/bash
#
# Fanuc-style G71.3 followed by G70.3 finishing on the same P-Q range.

rs274 -g g71.3-g70.3.ngc | awk '{$1=""; print}' > result &
pid=$!

count=10
while [ 0 -lt $count ] && kill -0 $pid > /dev/null 2>&1 ; do
    sleep 1
    count=$((count - 1))
done

if kill -0 $pid > /dev/null 2>&1; then
    kill -9 $pid
    echo "E: g71.3-g70.3.ngc seems stuck, killing"
    exit 1
fi

wait "$pid"
exit $?
