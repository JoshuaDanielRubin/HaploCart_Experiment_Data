#!/bin/bash

thousand_genomes_dir="../src/simulations/thousand_genomes"

#for id in NA19661 HG00473 HG01051 HG02666 HG03112 NA18510 NA19036 NA20518; do
#    for bam in $thousand_genomes_dir/$id/*.bam; do
#        echo $bam
#        echo -n $bam >> ../data/1kdepths.txt
#        echo -n " " >> ../data/1kdepths.txt
#        samtools view -h $bam chrM | samtools depth -a - | awk '{sum+=$3} END {print sum/16571}' >> ../data/1kdepths.txt
#    done
#done

for bam in ../src/simulations/alignments/*.bam; do
    echo $bam
    echo -n $bam >> ../data/fastq_no_numt_sim_depths.txt
    echo -n " " >>  ../data/fastq_no_numt_sim_depths.txt
        depth=$(samtools depth -a $bam | awk '{sum+=$3} END {print sum/16571}')
        echo $depth >> ../data/fastq_no_numt_sim_depths.txt
    done

#for bam in ../src/simulations/alignments/numt/*.bam; do
#    echo $bam
#    echo -n $bam >> ../data/fastq_with_numt_sim_depths.txt
#    echo -n " " >>  ../data/fastq_with_numt_sim_depths.txt
#        depth=$(samtools depth -a $bam | awk '{sum+=$3} END {print sum/16571}')
#        echo $depth >> ../data/fastq_with_numt_sim_depths.txt
#    done


