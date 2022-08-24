
get_count_per_bin.sh - Get counts per window of mean coverage depth for plotting purposes.

get_depths.sh - Get average depth of coverage of files for use in downsampling experiments.

hsd2fasta.py - Convert HSD to FASTA.

circularize_paths.py - Circularize all paths in a graph in ODGI format.

get_haplogrep_preds.sh - Obtain Haplogrep2 predictions using the RCRS-based Phylotree version 17.

mask_fasta.py - Mask a certain number of bases in a FASTA file.

create_pangenome_mapping.sh - Use ODGI to obtain a mapping of node ids to pangenome bases.

parse_pangenome_mapping.sh - Parse the output of create_pangenome_mapping.sh to obtain a TSV file.

plot_downsample.py - Create plots for downsampling experiments.

plot_mask.py - Plot distribution of edit distances for the masking experiment.

plot_mask_correctness.py - Plot proportion exactly correct in the masking experiment.

plot_qual.py - Plot reported confidences in predictions.

run_pggb.sh - Construct the graph.

score.py - Routines for scoring predictions in benchmarking experiments.

wrapper.py - Wrapper script for the web app.

write_hsd.py - Given a set of variants with respect to the RSRS, create HSD files (note this was not strictly
               necessary, we could have gone straight from variants to FASTA instead). 

benchmark.py - scripts to run benchmarking experiments.

run_tests.py - Run unit tests.
