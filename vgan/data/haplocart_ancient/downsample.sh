#!/bin/bash

# Input parameter
bam_dir=$1 # Path to directory containing BAM files

# Create directories for the subsampled BAM files at different depths
mkdir -p subsampled/1x
#mkdir -p subsampled/2x
#mkdir -p subsampled/0.5x

# Function to subsample a BAM file at a given depth
subsample_bam() {
  bam_file=$1
  chr=$2
  target_depth=$3

  # Calculate current depth of coverage on the mitochondria
  depth=$(samtools depth -r $chr $bam_file | awk '{sum+=$3} END {print sum/NR}')

  # Calculate the fraction of reads to subsample
  #fraction=$(echo "$target_depth/$depth" | bc -l)

  # Subsample the BAM file
  #filename=$(basename $bam_file)
  #out_dir="subsampled/${target_depth}x"
  #out_file="${filename%%.bam}_${target_depth}x.bam"
  #nice samtools view -bs $fraction $bam_file -@ 60 > $out_dir/$out_file

  # Calculate the depth of coverage of the subsampled BAM file
  #subsampled_depth=$(samtools depth $out_dir/$out_file | awk '{sum+=$3} END {print sum/NR}')

  echo "Depth of coverage of $bam_file on the mitochondria is $depth"
  #echo "Subsampled BAM file to a depth of $subsampled_depth"
  #echo "Subsampled BAM file saved as $out_dir/$out_file"
}

# Export the subsample_bam function so it can be used by GNU Parallel
export -f subsample_bam

# Loop through all BAM files in the input directory
for bam_file in ${bam_dir}/*.bam; do

  # Check if mitochondrial chromosome is called "chrM" or "MT"
  if samtools idxstats $bam_file | cut -f1 | grep -q "chrM"; then
    chr="chrM"
  elif samtools idxstats $bam_file | cut -f1 | grep -q "MT"; then
    chr="MT"
  else
    echo "Mitochondrial chromosome not found in $bam_file"
    continue
  fi

  # Subsample the BAM file using GNU Parallel
  filename=$(basename $bam_file)
  echo "Subsampling $filename..."

  # Subsample at target depth of 1x
  nice parallel -j 20 subsample_bam {} $chr 1 ::: $bam_file

  # Subsample at target depth of 2x
  #nice parallel -j 20 subsample_bam {} $chr 2 ::: $bam_file

  # Subsample at target depth of 0.5x
  #nice parallel -j 20 subsample_bam {} $chr 0.5 ::: $bam_file

done

