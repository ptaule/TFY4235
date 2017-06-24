#!/usr/bin/env bash

make all
binary=./percolation.prog

for i in `seq 5 10`;
do
    time ${binary} <<< $((100* ${i}))
done
