#!/bin/bash

# This script runs the entire data processing and import pipeline in the correct order.
# It includes prerequisite downloads and uses robust paths for executables.

# Exit immediately if a command exits with a non-zero status.
set -e

echo "--- Starting Full Data Import and Processing Pipeline ---"

# --- Step 1: Pre-process VCF Data ---

echo "[1/11] Filtering for significant variants..."
bcftools view -i 'FORMAT/LP > 7.3' data/ukb-d-2395_1.vcf.gz -o results/significant_hg37.vcf

echo "[2/11] Performing LiftOver from hg37 to hg38..."
echo "Downloading prerequisite files for LiftOver (if they don't exist)..."
wget -nc -P downloads/ [http://hgdownload.soe.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz](http://hgdownload.soe.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz)
wget -nc -P downloads/ [http://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz](http://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz)
gunzip -f downloads/hg19ToHg38.over.chain.gz || true
gunzip -f downloads/hg38.fa.gz || true

echo "Running CrossMap..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate gwas-env
$CONDA_PREFIX/bin/CrossMap vcf downloads/hg19ToHg38.over.chain results/significant_hg37.vcf downloads/hg38.fa results/gwas_hg38.vcf


echo "[3/11] Re-annotating with VEP..."
echo "PAUSED: Please annotate 'results/gwas_hg38.vcf' with VEP (GRCh38) and save as 'results/egIwc7NRt4hou5yo.txt'"
read -p "Press [Enter] to continue once the annotation file is ready..."

# --- Step 2: Build the Knowledge Graph ---
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

# --- Step 3: Enrich Graph with Unstructured Data ---
echo "[9/11] Fetching PubMed abstracts..."
python 6_pubmed_fetcher.py

echo "[10/11] Extracting entities from PubMed abstracts and enriching graph..."
python 7_ner_importer.py

echo "[11/11] Reconciling extracted entities with known graph nodes..."
python 8_ner_reconciliation.py

# --- Step 4: Build the Vector Database ---
echo "--- Building Vector Database ---"
echo "Creating embeddings and storing in ChromaDB..."
python 9_create_embeddings.py

echo "--- Pipeline Complete! ---"
echo "You can now run 'python query_engine.py' to chat with your knowledge graph."

