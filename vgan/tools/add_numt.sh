#!/bin/bash
IFS=
newline=$'\n'
fastqin=$1
numt=$2
rate=$3
out=$4

Nlines=$(zcat $fastqin | wc -l)
Nreads=$((Nlines/4))

ret="$(zcat $fastqin)"
for ((c=1; c<=Nreads; c++)); do
    if ((RANDOM % $rate == 0)); then
        random_numt=$((RANDOM % 250))
        numt_to_add=$(zcat $numt | tail -n $(($random_numt*4)) | head -n 4)
        ret+=$newline;
        ret+=$numt_to_add;
                             fi
                             done

echo $ret | gzip > $out
