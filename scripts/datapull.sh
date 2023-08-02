#!/bin/bash

# Set the input and output directory paths
output_dir="/home/jupyter-tettevie/advir/output"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Read the accession list file line by line
while IFS= read -r accession; do
 
 # Run sratools command to retrieve the data
  prefetch "$accession" -O "$output_dir"
  
  # Run fastq-dump command to convert SRA to FASTQ
  fastq-dump --outdir "$output_dir" --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 "$accession"

  # Perform quality control using FastQC
  fastqc "$output_dir/$accession.fastq" -o "$output_dir/fastqc_output"

  # Perform adapter trimming and quality filtering using trim_galore
  trim_galore --fastqc --illumina --length 36 --output_dir "$output_dir" "$output_dir/$accession.fastq"

done <<EOF
ERR9452447
SRR16303551
EOF