#!/bin/bash

# This script runs the entire data processing and import pipeline in the correct order.
# It now includes all data sources, including ENCODE regulatory elements.

# Exit immediately if a command exits with a non-zero status.
set -e

echo "--- Starting Full Data Import and Processing Pipeline ---"

# --- Step 0: Setup Directories ---
echo "[0/14] Creating project directories if they do not exist..."
mkdir -p data
mkdir -p models
mkdir -p results
mkdir -p downloads
mkdir -p neo4j/data

# --- Initialize Conda in the script's shell ---
# This ensures that tools installed in the active environment are found.
echo "Initializing Conda for this script session..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate gwas-env

# --- Step 1: Pre-process VCF Data ---

echo "[1/14] Filtering for significant variants..."
bcftools view -i 'FORMAT/LP > 7.3' data/ukb-d-2395_1.vcf.gz -o results/significant_hg37.vcf

echo "[2/14] Performing LiftOver from hg37 to hg38..."
wget -nc -P downloads/ http://hgdownload.soe.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
wget -nc -P downloads/ http://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
gunzip -f downloads/hg19ToHg38.over.chain.gz || true
gunzip -f downloads/hg38.fa.gz || true
echo "Running CrossMap..."
$CONDA_PREFIX/bin/CrossMap vcf downloads/hg19ToHg38.over.chain results/significant_hg37.vcf downloads/hg38.fa results/gwas_hg38.vcf


# --- Step 2: Annotate with VEP using Docker ---
echo "[3/14] Annotating with VEP via Docker container..."

# Permission checks for Docker
echo "Verifying write permissions for critical directories..."
sudo chmod -R 777 ~/.vep results neo4j/data
echo "Permissions set."


# This command runs VEP inside its official Docker container.
docker run --rm -v $(pwd)/results:/opt/vep/data -v ~/.vep:/opt/vep/.vep ensemblorg/ensembl-vep \
    vep -i /opt/vep/data/gwas_hg38.vcf -o /opt/vep/data/egIwc7NRt4hou5yo.txt \
    --cache \
    --dir /opt/vep/.vep \
    --assembly GRCh38 \
    --fork 12 \
    --force_overwrite \
    --tab \
    --distance 2000 \
    --symbol \
    --biotype \
    --uniprot \
    --fields "Uploaded_variation,Location,Allele,Consequence,IMPACT,SYMBOL,REF_ALLELE,Gene,GO,SWISSPROT,TREMBL,CHROM,POS" \
    --plugin GO


# --- Step 3: Knowledge Graph Construction ---
echo "[4/14] Building base graph (Mutations -> Genes)..."
python 1_neo4j_base_importer.py

echo "[5/14] Importing Gene Ontology (GO) to link genes to biological processes..."
python 2_go_importer.py

echo "[6/14] Importing Human Phenotype Ontology (HPO) to link genes to phenotypes..."
python 3_hpo_importer.py

echo "[7/14] Importing Reactome pathways to link genes to biological pathways..."
python 4_reactome_importer.py

echo "[8/14] Importing ClinVar data to link variants to clinical significance..."
python 5_clinvar_importer.py

echo "[9/14] Importing DGIdb drug-gene interaction data..."
python 10_dgidb_importer.py

echo "[10/14] Importing STRING-DB protein-protein interaction data..."
python 11_ppi_importer.py

echo "[11/14] Importing ENCODE regulatory element data and linking to variants..."
python 12_encode_importer.py

# --- Step 4: Unstructured Data Enrichment and Vector DB Construction ---
echo "[12/14] Fetching PubMed abstracts for contextually relevant genes..."
python 6_pubmed_fetcher.py

echo "[13/14] Extracting entities (NER) from abstracts and adding to the graph..."
python 7_ner_importer.py

echo "[14/14] Reconciling extracted entities with known graph nodes..."
python 8_ner_reconciliation.py

echo "--- Building Vector Database ---"
echo "Creating embeddings and storing in ChromaDB..."
python 9_create_embeddings.py

echo "--- Pipeline Complete! ---"
echo "You can now run 'python query_engine.py' to chat with your knowledge graph."
