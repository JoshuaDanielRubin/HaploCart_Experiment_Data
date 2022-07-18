#!/bin/bash -xv

fasta_path="all.fa"
pggb_path="/home/ctools/pggb/pggb"

nice $pggb_path -i $fasta_path -t 63 -o synthetic -p 97 -s 5000 -H 5180 -n 5180 -k 300 -G 3079,3559


