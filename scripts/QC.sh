#!/bin/bash

# Specify the directory containing the FASTQ files
fastq_dir="/home/jupyter-tettevie/advir/output"

# Specify the directory to store the FastQC output
output_dir="/home/jupyter-tettevie/advir/output/fastqc_output"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through each FASTQ file in the directory
for fastq_file in "$fastq_dir"/*.fastq; do
    # Run FastQC on the current file
    fastqc "$fastq_file" --outdir "$output_dir"
done