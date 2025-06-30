# Bio-KG RAG: A Knowledge-Graph-Powered RAG System for GWAS Data

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/release/python-390/)

## Abstract

Genome-Wide Association Studies (GWAS) are instrumental in identifying genetic variants associated with complex traits and diseases. However, a significant challenge lies in interpreting these findings to understand the underlying biological mechanisms. This project presents a complete, reproducible pipeline for constructing a multi-layered biomedical knowledge graph from raw GWAS summary statistics and leveraging it within a sophisticated Retrieval-Augmented Generation (RAG) system. By integrating structured ontological data (Gene Ontology, HPO, Reactome), clinical significance (ClinVar), protein-protein interactions (STRING-DB), regulatory elements (ENCODE), and drug-gene interactions (DGIdb) with unstructured knowledge extracted from scientific literature (PubMed), this system allows researchers to ask complex, natural language questions and receive synthesized, evidence-based answers, thereby bridging the gap between statistical genetic findings and actionable biological insight.

---

## Dataset

This pipeline was developed and tested using publicly available GWAS summary statistics from the IEU OpenGWAS project.

-   **Trait:** Hair/balding pattern: Pattern 1
-   **IEU OpenGWAS ID:** `ukb-d-2395_1`
-   **Source VCF File:** `ukb-d-2395_1.vcf.gz` (hg37/GRCh37 assembly)
-   **Link:** [https://gwas.mrcieu.ac.uk/datasets/ukb-d-2395_1/](https://gwas.mrcieu.ac.uk/datasets/ukb-d-2395_1/)

---

## System Architecture

This project implements a GraphRAG architecture composed of three core components:

1.  **Knowledge Graph (Neo4j):** A graph database that stores structured, interconnected data. This includes GWAS variants, their affected genes, protein-protein interactions, biological pathways, associated phenotypes, known drug targets, regulatory elements, and entities extracted from literature.
2.  **Vector Database (ChromaDB):** A vector store containing embeddings of scientific abstracts from PubMed, enabling fast semantic search.
3.  **Query Engine (Local LLM):** A `Phi-3-medium` language model orchestrated by a Python script that uses a biomedical NER model to understand queries, retrieves context from both the knowledge graph and the vector database, and synthesizes this information to generate a comprehensive, evidence-based answer.

---

## Project Structure

The project is organized into a modular pipeline, with each script performing a distinct data ingestion or processing task.

```
/
├── data/
│   └── ukb-d-2395_1.vcf.gz         # Original GWAS VCF file (hg37)
├── models/
│   └── Phi-3-medium-128k-instruct-Q6_K_L.gguf  # The LLM file
├── results/                        # Directory for generated output files
├── downloads/                      # Directory for downloaded prerequisite files
├── neo4j/
│   └── data/                       # Neo4j database files (mounted volume)
│
├── docker-compose.yml              # Neo4j service configuration
├── README.md                       # This documentation file
├── requirements.txt                # A list of all required Python packages
|
├── run_pipeline.sh                 # Master script to run the full data pipeline
├── clean.sh                        # Utility script to clean the directory
|
├── 1_neo4j_base_importer.py        # Loads base mutations and genes
├── 1.5_gwas_context_importer.py    # Adds GWAS study metadata to the graph
├── 2_go_importer.py                # Loads Gene Ontology data and relationships
├── 3_hpo_importer.py               # Loads Human Phenotype Ontology data
├── 4_reactome_importer.py          # Loads Reactome pathways and gene links
├── 5_clinvar_importer.py           # Loads ClinVar clinical significance data
├── 6_pubmed_fetcher.py             # Downloads abstracts from PubMed for relevant genes
├── 7_ner_importer.py               # Extracts entities from abstracts and adds to graph
├── 8_ner_reconciliation.py         # Reconciles extracted entities with known nodes
├── 9_create_embeddings.py          # Builds the vector database from abstracts
├── 10_dgidb_importer.py            # Loads drug-gene interaction data from DGIdb
├── 11_ppi_importer.py              # Loads protein-protein interactions from STRING-DB
├── 12_encode_importer.py           # Loads regulatory elements and cell-type data from ENCODE
|
└── query_engine.py                 # The final, interactive RAG query engine
```

---

## How to Run the System

### Step 0: Initial Setup

1.  **Conda Environment:** Create and activate a single conda environment for the project.
    ```bash
    conda create -n gwas-env python=3.9
    conda activate gwas-env
    ```

2.  **Install Dependencies:** Install all required packages.
    ```bash
    # Install command-line bioinformatics tools
    conda install -c bioconda bcftools crossmap bgzip tabix

    # Install Python packages
    pip install -r requirements.txt
    pip install scispacy==0.5.1 [https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.1/en_core_sci_lg-0.5.1.tar.gz](https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.1/en_core_sci_lg-0.5.1.tar.gz)
    
    # This final command requires system build tools like build-essential and cmake
    CMAKE_ARGS="-DGGML_CUDA=on" FORCE_CMAKE=1 pip install llama-cpp-python --force-reinstall --upgrade --no-cache-dir
    ```

3.  **Manual Data Downloads:** Before running the pipeline, you must manually download two files:
    * **GWAS VCF File:** Download the source `ukb-d-2395_1.vcf.gz` file and place it in the `data/` directory.
        * **Link:** [https://gwas.mrcieu.ac.uk/datasets/ukb-d-2395_1/](https://gwas.mrcieu.ac.uk/datasets/ukb-d-2395_1/)
    * **ENCODE Cell-Type Specificity Data:** Download the DNase-Z-score matrix from the SCREEN portal:
        * **URL:** [https://screen.encodeproject.org/downloads](https://screen.encodeproject.org/downloads)
        * **File:** Find the "Human (GRCh38/hg38)" section and download the file listed as "DNase Z-score matrix". It will have a filename like `ENCFF833CSA.txt.gz`.
        * Place this file directly into the `downloads/` directory.

4.  **Set up VEP Docker Environment & Permissions:**
    Pull the latest VEP Docker image and set the necessary directory permissions.
    ```bash
    docker pull ensemblorg/ensembl-vep:latest
    mkdir -p ~/.vep results neo4j/data
    sudo chmod -R 777 ~/.vep results neo4j/data
    ```

5.  **Clean the Directory (Optional):** To start completely fresh, run the `clean.sh` script.
    ```bash
    bash clean.sh
    ```

### Step 1: Start the Database Service

Start the Neo4j Docker container in the background.
```bash
docker-compose up -d neo4j
```

### Step 2: Run the Full Data Pipeline

Execute the master script. This will run the entire, fully automated data processing and annotation pipeline.
```bash
bash run_pipeline.sh
```
*Note: The script now automatically handles all data downloads and processing steps, aside from the two manual downloads listed above.*

### Step 3: Query Your Knowledge Graph

Once the pipeline is complete, start the interactive RAG engine.
```bash
python query_engine.py
```
