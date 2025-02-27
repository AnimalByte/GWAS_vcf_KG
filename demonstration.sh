#!/usr/bin/env bash

# Exit immediately if a command exits with a non-zero status.
set -e

# Enable debug mode (optional, uncomment for troubleshooting)
# set -x

# Define directories
DATA_DIR="data"
RESULTS_DIR="results"
ANNOTATIONS_DIR="annotations"

# Create necessary directories if they don't exist
mkdir -p "$DATA_DIR" "$RESULTS_DIR" "$ANNOTATIONS_DIR"

# Step 1: Extract GWAS Summary Statistics from the VCF
echo "Step 1: Extracting GWAS Summary Statistics from the VCF..."

# Write header to gwas_summary_stats.txt
echo -e "CHROM\tPOS\tID\tREF\tALT\tAF\tES\tSE\tLP" > "${DATA_DIR}/gwas_summary_stats.txt"

# Corrected bcftools query command
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\t%FORMAT/ES\t%FORMAT/SE\t%FORMAT/LP\n' "${DATA_DIR}/ukb-b-12064.vcf.gz" >> "${DATA_DIR}/gwas_summary_stats.txt"

echo "Step 1 completed: GWAS summary statistics extracted to ${DATA_DIR}/gwas_summary_stats.txt"

# Step 2: Filter for Significant Variants (LP > -log10(5e-8))
echo "Step 2: Filtering for significant variants (LP > -log10(5e-8))..."

# Calculate the threshold using awk
threshold=$(awk 'BEGIN { print -log(5e-8)/log(10) }')
echo "Calculated threshold (LP > $threshold)"

# Filter significant variants
awk -v thresh="$threshold" 'NR > 1 && $9 >= thresh' "${DATA_DIR}/gwas_summary_stats.txt" > "${RESULTS_DIR}/significant_variants.txt"

echo "Step 2 completed: Significant variants saved to ${RESULTS_DIR}/significant_variants.txt"

# Step 3: Convert Significant Variants to VCF Format
echo "Step 3: Converting significant variants to VCF format..."

awk 'BEGIN {
  # Define VCF header
  print "##fileformat=VCFv4.2";
  print "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">";
  print "##INFO=<ID=LP,Number=1,Type=Float,Description=\"-log10(p-value)\">";
  print "##FILTER=<ID=PASS,Description=\"All filters passed\">";
  print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
}
NR > 1 {  # Skip the header line
  print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t.\tPASS\tAF="$6";LP="$9
}' "${RESULTS_DIR}/significant_variants.txt" > "${RESULTS_DIR}/significant_variants.vcf"

echo "Step 3 completed: VCF file created at ${RESULTS_DIR}/significant_variants.vcf"

# Step 4: Annotate Significant Variants with VEP
echo "Step 4: Annotating significant variants with VEP..."

vep \
  --force_overwrite \
  --offline \
  --cache \
  --species homo_sapiens \
  --assembly GRCh37 \
  -i "${RESULTS_DIR}/significant_variants.vcf" \
  -o "${RESULTS_DIR}/vep_annotated_offline.txt" \
  --tab \
  --fields "Uploaded_variation,Location,Allele,Gene,Feature,Feature_type,Consequence,IMPACT,GO_terms,GO_biological_process,GO_cellular_component,GO_molecular_function" \
  --symbol \
  --canonical \
  --plugin GO \
  --fork 6



echo "Step 4 completed: VEP annotation saved to ${RESULTS_DIR}/vep_annotated.txt"

# Step 5: Map Variants to Genes
echo "Step 5: Mapping variants to genes..."

# Extract SNP, Location, and Gene columns (assuming tab-separated)
awk 'BEGIN {FS="\t"; OFS="\t"} NR > 1 { print $1, $2, $3 }' "${RESULTS_DIR}/vep_annotated.txt" > "${ANNOTATIONS_DIR}/snp_gene_map.txt"

echo "Step 5 completed: SNP to gene mapping saved to ${ANNOTATIONS_DIR}/snp_gene_map.txt"

# Confirm pipeline completion
echo "Pipeline completed successfully. Outputs are available in the '${RESULTS_DIR}/' and '${ANNOTATIONS_DIR}/' directories."
