#!/bin/bash

# Define input and output folder paths
input_folder="data/vcf_input"
output_folder="data/vcf_annotated"

# Create output folder if it does not exist
mkdir -p "$output_folder"

# Loop through all .freebayes.vcf files in the input folder
for vcf_file in "$input_folder"/*.freebayes.vcf; do
    if [ ! -e "$vcf_file" ]; then
        echo "No .freebayes.vcf files found in input folder!"
        break
    fi

    # Extract file name without path
    filename=$(basename "$vcf_file")

    # Define output file path
    output_file="$output_folder/${filename%.vcf}.annotated.vcf"

    # Run VEP annotation
    tools/ensembl-vep/vep \
        --input_file "$vcf_file" \
        --output_file "$output_file" \
        --format vcf \
        --vcf \
        --species homo_sapiens \
        --assembly GRCh37 \
        --everything \
        --force_overwrite \
        --fasta reference/human_g1k_v37.fasta \
        --offline \
        --dir_cache tools/vep_cache

    echo "Annotation completed: $filename -> $output_file"
done

echo "All files have been annotated."

