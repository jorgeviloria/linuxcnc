#!/bin/bash
#
# G71.1 + G71.2: the envelope connected to Z0 disappears before the groove
# bottom, and G71.2 must still cut the pocket.

rs274 -g g71.2-pocket-only.ngc | awk '{$1=""; print}' > result &
pid=$!

count=10
while [ 0 -lt $count ] && kill -0 $pid > /dev/null 2>&1 ; do
    sleep 1
    count=$((count - 1))
done

if kill -0 $pid > /dev/null 2>&1; then
    kill -9 $pid
    echo "E: g71.2-pocket-only.ngc seems stuck, killing"
    exit 1
fi

wait "$pid"
exit $?
