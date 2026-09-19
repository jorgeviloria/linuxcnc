#!/bin/bash
#
# Fanuc-style G71.3: the scan must reach the real bottom of an arc whose
# interior X is lower than its end points.

rs274 -g g71.3-arc-bottom.ngc | awk '{$1=""; print}' > result &
pid=$!

count=10
while [ 0 -lt $count ] && kill -0 $pid > /dev/null 2>&1 ; do
    sleep 1
    count=$((count - 1))
done

if kill -0 $pid > /dev/null 2>&1; then
    kill -9 $pid
    echo "E: g71.3-arc-bottom.ngc seems stuck, killing"
    exit 1
fi

wait "$pid"
exit $?
