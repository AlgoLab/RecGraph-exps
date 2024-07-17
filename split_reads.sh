#!/bin/bash

# Define the genes and error levels
GENES=$(cat genes_HLA_short.txt)
ERRORS=("0" "3" "5" "10" "rec")

# Iterate over each gene and error level
for GENE in $GENES; do
    for ERR in "${ERRORS[@]}"; do
        INPUT_FILE="output/HLA/genes/${GENE}/reads_${ERR}.fa"
        OUTPUT_DIR="output/HLA/genes/${GENE}/reads_${ERR}_split"
        
        # Create the output directory if it doesn't exist
        mkdir -p $OUTPUT_DIR
        
        # Split the input file into individual read files
        awk -v OUT_DIR="$OUTPUT_DIR" '
        /^>/ {
            if (out_file) close(out_file)
            read_count++
            out_file = sprintf("%s/read_%d.fa", OUT_DIR, read_count)
        }
        { print > out_file }
        ' $INPUT_FILE
        rm -f $INPUT_FILE

    done
done
