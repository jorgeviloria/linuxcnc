#!/bin/bash
#
# Two G71.3/G70.3 pairs in one file: the first G70.3 falls back to a
# search (the cache holds the second range), the second uses the cache.

rs274 -g g71.3-two-ranges.ngc | awk '{$1=""; print}' > result &
pid=$!

count=10
while [ 0 -lt $count ] && kill -0 $pid > /dev/null 2>&1 ; do
    sleep 1
    count=$((count - 1))
done

if kill -0 $pid > /dev/null 2>&1; then
    kill -9 $pid
    echo "E: g71.3-two-ranges.ngc seems stuck, killing"
    exit 1
fi

wait "$pid"
exit $?
