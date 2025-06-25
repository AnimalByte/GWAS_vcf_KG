# Bio-KG RAG: A Knowledge-Graph-Powered RAG System for GWAS Data

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/release/python-390/)

## Abstract

Genome-Wide Association Studies (GWAS) are instrumental in identifying genetic variants associated with complex traits and diseases. However, a significant challenge lies in interpreting these findings to understand the underlying biological mechanisms. This project presents a complete, reproducible pipeline for constructing a multi-layered biomedical knowledge graph from raw GWAS summary statistics and leveraging it within a sophisticated Retrieval-Augmented Generation (RAG) system. By integrating structured ontological data (GO, HPO, Reactome, ClinVar) with unstructured knowledge extracted from scientific literature (PubMed), this system allows researchers to ask complex, natural language questions and receive synthesized, evidence-based answers, thereby bridging the gap between statistical genetic findings and actionable biological insight.

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

## Methodology: The Data Pipeline

The pipeline is executed by the `run_pipeline.sh` script, which orchestrates the following steps in order:

#### Pre-processing
1.  **Filtering for Significance:** The raw GWAS VCF file is first filtered using `bcftools` to retain only genome-wide significant variants (P < 5x10⁻⁸).
2.  **Genomic Coordinate Harmonization:** The coordinates of the significant variants are converted ("lifted over") from the hg37/GRCh37 assembly to the modern hg38/GRCh38 assembly using `CrossMap`.
3.  **Functional Annotation:** The harmonized VCF file is annotated using the Ensembl Variant Effect Predictor (VEP) to identify affected genes and predict functional consequences.

#### Knowledge Graph Construction
The core of the system is a Neo4j graph database, populated by a series of modular Python scripts:

4.  **Base Graph Import (`1_neo4j_base_importer.py`):** Creates `:Mutation` and `:Gene` nodes.
5.  **Ontology Integration (`2_go_importer.py`, `3_hpo_importer.py`):** Imports the Gene Ontology (GO) and Human Phenotype Ontology (HPO), creating `:GO_Term` and `:Phenotype` nodes and their rich hierarchies.
6.  **Pathway Integration (`4_reactome_importer.py`):** Adds `:Pathway` nodes from the Reactome database.
7.  **Clinical Significance (`5_clinvar_importer.py`):** Enriches `:Mutation` nodes with clinical significance data from ClinVar.

#### Unstructured Data Integration
8.  **Literature Retrieval (`6_pubmed_fetcher.py`):** Fetches relevant scientific abstracts from PubMed for genes in the graph.
9.  **Knowledge Extraction (`7_ner_importer.py`):** Uses a GPU-accelerated `scispaCy` NER model to identify biomedical entities in abstracts and adds them to the graph.
10. **Entity Reconciliation (`8_ner_reconciliation.py`):** Links the unstructured and structured worlds by merging entities found in literature with the official ontology nodes in the graph.

#### RAG System Preparation
11. **Vector Embedding (`9_create_embeddings.py`):** Uses a PubMedBERT model to create vector embeddings for all abstracts and stores them in a local ChromaDB database.

---

## How to Run the System

### Step 0: Initial Setup

1.  **Conda Environment:** Create and activate the `gwas-env` conda environment.
    ```bash
    conda create -n gwas-env python=3.9
    conda activate gwas-env
    ```

2.  **Dependencies:** Install all required Python packages using the provided `requirements.txt` file. It's recommended to install `llama-cpp-python` with GPU support as a separate, final step.
    ```bash
    pip install -r requirements.txt
    CMAKE_ARGS="-DGGML_CUDA=on" FORCE_CMAKE=1 pip install llama-cpp-python --force-reinstall --upgrade --no-cache-dir
    ```
    *(Note: You will also need system build tools like `build-essential` and `cmake` installed).*

3.  **Clean the Directory (Optional):** To start completely fresh, run the `clean.sh` script.
    ```bash
    bash clean.sh
    ```

### Step 1: Start the Database Service

Start the Neo4j Docker container in the background.
```bash
docker-compose up -d neo4j
```

### Step 2: Run the Full Data Pipeline

Execute the master script. This will run all data processing and import steps. **Note:** The script will pause for manual VEP annotation.
```bash
bash run_pipeline.sh
```

### Step 3: Query Your Knowledge Graph

Once the pipeline is complete, start the interactive RAG engine to ask questions.
```bash
python query_engine.py
```
