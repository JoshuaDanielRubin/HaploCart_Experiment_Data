#!/bin/bash

# Replace /path/to/bam/files with the path to your bam files directory
BAM_DIR=.

# Create a new directory to hold the individual bam files and their indexes
mkdir -p ${BAM_DIR}/individual_bams

# Loop through each bam file in the bam directory
for bam_file in ${BAM_DIR}/*.bam; do

    # Extract the base filename (without extension) of the bam file
    base_filename=$(basename ${bam_file} .bam)

    # Create a new directory for this bam file and its index
    mkdir -p ${BAM_DIR}/individual_bams/${base_filename}

    # Move the bam file and its corresponding index into the new directory
    cp ${bam_file} ${BAM_DIR}/individual_bams/${base_filename}/
    cp ${bam_file}.bai ${BAM_DIR}/individual_bams/${base_filename}/

done

