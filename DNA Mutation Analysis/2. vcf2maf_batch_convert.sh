#!/bin/bash

# Set input and output folders
input_folder="input_vcf_annotated"  # Folder containing VEP-annotated VCFs
output_folder="output_maf_files"    # Output folder for MAF files

# Create output directory if it doesn't exist
mkdir -p "$output_folder"

# Loop over all VCF files in the input folder
for vcf_file in "$input_folder"/*.vcf; do
    # Extract base filename (without path or .vcf extension)
    filename=$(basename "$vcf_file" .vcf)
    
    # Define output MAF file path
    output_maf="$output_folder/${filename}.maf"

    # Convert using vcf2maf
    perl <path_to_vcf2maf>/vcf2maf.pl \
        --input-vcf "$vcf_file" \
        --output-maf "$output_maf" \
        --tumor-id "$filename" \
        --ref-fasta <path_to_reference>/human_g1k_v37.fasta \
        --vep-path <path_to_vep> \
        --vep-data <path_to_vep_cache> \
        --tmp-dir "$output_folder/tmp"

    echo "Converted $vcf_file -> $output_maf"
done

echo "All VCFs converted to MAF!"
