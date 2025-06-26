# Bio-KG RAG: A Knowledge-Graph-Powered RAG System for GWAS Data

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/release/python-390/)

## Abstract

Genome-Wide Association Studies (GWAS) are instrumental in identifying genetic variants associated with complex traits and diseases. However, a significant challenge lies in interpreting these findings to understand the underlying biological mechanisms. This project presents a complete, reproducible pipeline for constructing a multi-layered biomedical knowledge graph from raw GWAS summary statistics and leveraging it within a sophisticated Retrieval-Augmented Generation (RAG) system. By integrating structured ontological data (GO, HPO, Reactome, ClinVar) with unstructured knowledge extracted from scientific literature (PubMed), this system allows researchers to ask complex, natural language questions and receive synthesized, evidence-based answers, thereby bridging the gap between statistical genetic findings and actionable biological insight.

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

1.  **Knowledge Graph (Neo4j):** A graph database that stores structured, interconnected data. This includes GWAS variants, their affected genes, biological pathways, associated phenotypes from clinical ontologies, and entities extracted from literature. This structured backbone allows for precise, multi-hop queries to uncover complex relationships.

2.  **Vector Database (ChromaDB):** A vector store containing embeddings of scientific abstracts from PubMed. This component enables fast and efficient **semantic search**, allowing the system to find documents based on conceptual meaning rather than just keywords. Embeddings are generated using a domain-specific PubMedBERT model.

3.  **Query Engine (Local LLM):** A large language model (`Phi-3-medium`) orchestrated by a Python script. It acts as the "brain" of the system, using a Named Entity Recognition (NER) model to understand user queries, retrieving context from both the knowledge graph and the vector database, and synthesizing this information to generate a comprehensive, evidence-based answer.

---

## Project Structure

The project is organized into a modular pipeline, with each script performing a distinct data ingestion or processing task.

```
/
├── data/
│   └── ukb-d-2395_1.vcf.gz         # Original GWAS VCF file (hg37)
├── models/
│   └── Phi-3-medium-128k-instruct-Q5_K_M.gguf  # The LLM file
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
├── 2_go_importer.py                # Loads Gene Ontology data and relationships
├── 3_hpo_importer.py               # Loads Human Phenotype Ontology data
├── 4_reactome_importer.py          # Loads Reactome pathways and gene links
├── 5_clinvar_importer.py           # Loads ClinVar clinical significance data
├── 6_pubmed_fetcher.py             # Downloads abstracts from PubMed for relevant genes
├── 7_ner_importer.py               # Extracts entities from abstracts and adds to graph
├── 8_ner_reconciliation.py         # Reconciles extracted entities with known nodes
├── 9_create_embeddings.py          # Builds the vector database from abstracts
|
└── query_engine.py                 # The final, interactive RAG query engine
```

---

## How to Run the System

### Step 0: Initial Setup

1.  **Conda Environment:** Create a single conda environment for the project.
    ```bash
    conda create -n gwas-env python=3.9
    conda activate gwas-env
    ```

2.  **Install Dependencies:** Install all required packages.
    ```bash
    # Install command-line bioinformatics tools
    conda install -c bioconda bcftools crossmap

    # Install Python packages
    pip install -r requirements.txt
    pip install scispacy==0.5.1 [https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.1/en_core_sci_lg-0.5.1.tar.gz](https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.1/en_core_sci_lg-0.5.1.tar.gz)
    
    # This final command requires system build tools like build-essential and cmake
    CMAKE_ARGS="-DGGML_CUDA=on" FORCE_CMAKE=1 pip install llama-cpp-python --force-reinstall --upgrade --no-cache-dir
    ```

3.  **Set up VEP Docker Environment (One-time only):**
    This involves pulling the VEP Docker image and creating a local directory to store its cache files.
    ```bash
    # Pull the latest official VEP Docker image
    docker pull ensemblorg/ensembl-vep:latest

    # Create a local directory for VEP data
    mkdir -p ~/.vep

    # Run the VEP installer inside the container to download cache and plugin data
    # This mounts your local .vep directory into the container for persistent storage.
    docker run -t -i -v ~/.vep:/opt/vep/.vep ensemblorg/ensembl-vep \
      INSTALL.pl -a cfp -s homo_sapiens -y GRCh38 --PLUGINS Geno2MP,ClinPred
    ```

4.  **Clean the Directory (Optional):** To start completely fresh, run the `clean.sh` script.
    ```bash
    bash clean.sh
    ```

### Step 1: Start the Database Service

Start the Neo4j Docker container in the background.
```bash
docker-compose up -d neo4j
```

### Step 2: Run the Full Data Pipeline

Execute the master script. This will run the entire, fully automated data processing pipeline.
```bash
bash run_pipeline.sh
```

### Step 3: Query Your Knowledge Graph

Once the pipeline is complete, activate your environment and start the interactive RAG engine.
```bash
conda activate gwas-env
python query_engine.py
```
