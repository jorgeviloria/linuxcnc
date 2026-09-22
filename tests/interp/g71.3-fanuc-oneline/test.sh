#!/bin/bash
#
# Fanuc-style G71.3 Type I, one-line form (I = depth of cut).

rs274 -g g71.3-oneline.ngc | awk '{$1=""; print}' > result &
pid=$!

count=10
while [ 0 -lt $count ] && kill -0 $pid > /dev/null 2>&1 ; do
    sleep 1
    count=$((count - 1))
done

if kill -0 $pid > /dev/null 2>&1; then
    kill -9 $pid
    echo "E: g71.3-oneline.ngc seems stuck, killing"
    exit 1
fi

wait "$pid"
exit $?
