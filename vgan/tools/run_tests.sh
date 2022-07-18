#!/bin/bash

rm -f ../test/output_files/*.txt;
for test in plus_strand_perfect_match plus_strand_mismatch plus_strand_insertion_in_read plus_strand_deletion_in_read plus_strand_softclip \
            minus_strand_perfect_match minus_strand_mismatch minus_strand_insertion_in_read minus_strand_deletion_in_read minus_strand_softclip \
            interleaved multifasta multifasta_zipped check_thread_too_many check_thread_zero check_thread_minus_one load check_nodevec check_graph rmdup gam \
            fq_single_rcrs consensus_rcrs zipped_consensus zipped_paired_fastq tough_consensus another_consensus custom_posterior_output \
            missing_input_consensus missing_input_fq1 missing_input_fq2 missing_input_gam invalid_bep1 invalid_bep2 valid_bep; do
    echo "Running test: " $test;
    ./../bin/test --run_test=$test;
done
