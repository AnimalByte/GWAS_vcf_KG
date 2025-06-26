#!/bin/bash

# This script runs the entire data processing and import pipeline in the correct order.
# It now uses the VEP Docker container for annotation, making the process more robust.

# Exit immediately if a command exits with a non-zero status.
set -e

echo "--- Starting Full Data Import and Processing Pipeline ---"

# --- Initialize Conda in the script's shell ---
# This ensures that tools installed in the active environment are found.
echo "Initializing Conda for this script session..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate gwas-env

# --- Step 1: Pre-process VCF Data ---

echo "[1/11] Filtering for significant variants..."
bcftools view -i 'FORMAT/LP > 7.3' data/ukb-d-2395_1.vcf.gz -o results/significant_hg37.vcf

echo "[2/11] Performing LiftOver from hg37 to hg38..."
# Download prerequisite files for LiftOver if they don't exist
echo "Downloading prerequisite files for LiftOver (if they don't exist)..."
wget -nc -P downloads/ http://hgdownload.soe.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
wget -nc -P downloads/ http://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
# Ensure the files are unzipped. The '|| true' prevents errors if the file is already unzipped.
gunzip -f downloads/hg19ToHg38.over.chain.gz || true
gunzip -f downloads/hg38.fa.gz || true

# Run CrossMap using the direct path from the conda environment
echo "Running CrossMap..."
$CONDA_PREFIX/bin/CrossMap vcf downloads/hg19ToHg38.over.chain results/significant_hg37.vcf downloads/hg38.fa results/gwas_hg38.vcf


# --- Step 2: Annotate with VEP using Docker ---
echo "[3/11] Annotating with VEP via Docker container..."

# *** NEW: Add permission checks for both directories needed by Docker ***
echo "Verifying write permissions for the '~/.vep' cache directory..."
if ! touch "$HOME/.vep/permission_test" 2>/dev/null; then
    echo "------------------------------------------------------------------"
    echo "ERROR: Permission denied."
    echo "The script cannot write to the VEP cache directory ('~/.vep')."
    echo "Please run the following command to fix the permissions, then"
    echo "re-run this script:"
    echo
    echo "sudo chmod -R 777 ~/.vep"
    echo "------------------------------------------------------------------"
    exit 1
fi
rm "$HOME/.vep/permission_test"
echo "VEP cache permissions are OK."


echo "Verifying write permissions for the 'results' directory..."
if ! touch "results/permission_test" 2>/dev/null; then
    echo "------------------------------------------------------------------"
    echo "ERROR: Permission denied."
    echo "The script cannot write to the './results' directory."
    echo "Please run the following command to fix the permissions, then"
    echo "re-run this script:"
    echo
    echo "sudo chmod -R 777 results"
    echo "------------------------------------------------------------------"
    exit 1
fi
rm "results/permission_test"
echo "Results directory permissions are OK."


# This command runs VEP inside its official Docker container.
# -v mounts the local 'results' directory into the container's working directory.
# -v also mounts the local VEP cache to prevent re-downloading.
docker run --rm -v $(pwd)/results:/opt/vep/data -v ~/.vep:/opt/vep/.vep ensemblorg/ensembl-vep \
    vep -i /opt/vep/data/gwas_hg38.vcf -o /opt/vep/data/egIwc7NRt4hou5yo.txt \
    --cache --offline --dir /opt/vep/.vep \
    --assembly GRCh38 \
    --fork 12 \
    --force_overwrite \
    --tab \
    --fields "Uploaded_variation,Location,Allele,Consequence,IMPACT,SYMBOL,REF_ALLELE,GO,SWISSPROT,TREMBL" \
    --plugin Geno2MP \
    --plugin ClinPred

# --- Step 3: Build the Knowledge Graph ---
echo "[4/11] Building base graph (Mutations and Genes)..."
python 1_neo4j_base_importer.py

echo "[5/11] Importing Gene Ontology (GO) data..."
python 2_go_importer.py

echo "[6/11] Importing Human Phenotype Ontology (HPO) data..."
python 3_hpo_importer.py

echo "[7/11] Importing Reactome Pathway data..."
python 4_reactome_importer.py

echo "[8/11] Importing ClinVar clinical significance data..."
python 5_clinvar_importer.py

# --- Step 4: Enrich Graph and Build Vector DB ---
echo "[9/11] Fetching PubMed abstracts..."
python 6_pubmed_fetcher.py

echo "[10/11] Extracting entities from abstracts and enriching graph..."
python 7_ner_importer.py

echo "[11/11] Reconciling extracted entities with known graph nodes..."
python 8_ner_reconciliation.py

echo "--- Building Vector Database ---"
echo "Creating embeddings and storing in ChromaDB..."
python 9_create_embeddings.py

echo "--- Pipeline Complete! ---"
echo "You can now run 'python query_engine.py' to chat with your knowledge graph."

