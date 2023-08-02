#!/bin/bash

# Set input and output directories
input_dir="/home/jupyter-tettevie/advir/output"
output_dir="/home/jupyter-tettevie/advir/output/trim_output"

# Iterate over all fastq files in the input directory
for file in "$input_dir"/*.fastq.gz; do
    # Get the base name of the file
    base_name=$(basename "$file" .fastq)
    
    # Set the output file name
    output_file="$output_dir/${base_name}_trimmed.fastq"
    
    # Run trim_galore with desired options
    trim_galore --quality 20 --length 50 --output_dir "$output_dir" "$file"
done