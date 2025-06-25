# Bio-KG RAG: A Knowledge-Graph-Powered RAG System for GWAS Data

## Abstract

Genome-Wide Association Studies (GWAS) are instrumental in identifying genetic variants associated with complex traits and diseases. However, a significant challenge lies in interpreting these findings to understand the underlying biological mechanisms. This project presents a complete, reproducible pipeline for constructing a multi-layered biomedical knowledge graph from raw GWAS summary statistics and leveraging it within a sophisticated Retrieval-Augmented Generation (RAG) system. By integrating structured ontological data (GO, HPO, Reactome, ClinVar) with unstructured knowledge extracted from scientific literature (PubMed), this system allows researchers to ask complex, natural language questions and receive synthesized, evidence-based answers, thereby bridging the gap between statistical genetic findings and actionable biological insight.

---

## Project Structure

The project is organized into a modular pipeline, with each script performing a distinct data ingestion or processing task.


/
|-- data/
|   |-- ukb-d-2395_1.vcf.gz         # Original GWAS VCF file (hg37)
|-- models/
|   |-- Phi-3-medium-128k-instruct-Q5_K_M.gguf  # The LLM file
|-- results/                        # All output files are generated here
|-- downloads/                      # Prerequisite files are downloaded here
|-- docker-compose.yml              # Neo4j service configuration
|
|-- 1_neo4j_base_importer.py        # Loads base mutations and genes
|-- 2_go_importer.py                # Loads Gene Ontology data and relationships
|-- 3_hpo_importer.py               # Loads Human Phenotype Ontology data
|-- 4_reactome_importer.py          # Loads Reactome pathways and gene links
|-- 5_clinvar_importer.py           # Loads ClinVar clinical significance data
|-- 6_pubmed_fetcher.py             # Downloads abstracts from PubMed for relevant genes
|-- 7_ner_importer.py               # Extracts entities from abstracts and adds to graph
|-- 8_ner_reconciliation.py         # Reconciles extracted entities with known nodes
|-- 9_create_embeddings.py          # Builds the vector database from abstracts
|
|-- query_engine.py                 # The final, interactive RAG query engine
|-- clean.sh                        # Utility script to clean the directory
|-- run_pipeline.sh                 # Master script to run the full pipeline
|-- README.md                       # This file


---

## Methodology: The Data Pipeline

The pipeline is executed by the `run_pipeline.sh` script, which orchestrates the following steps in order:

#### Pre-processing
1.  **Filtering for Significance:** The raw GWAS VCF file is first filtered using `bcftools` to retain only genome-wide significant variants (P < 5x10⁻⁸). This focuses the analysis on the most statistically relevant data.
2.  **Genomic Coordinate Harmonization:** The coordinates of the significant variants are converted ("lifted over") from the hg37/GRCh37 assembly to the modern hg38/GRCh38 assembly using `CrossMap`. This is a critical step to ensure compatibility with modern annotation databases.
3.  **Functional Annotation:** The harmonized VCF file is annotated using the Ensembl Variant Effect Predictor (VEP) to identify the affected genes, predict the functional consequence of each variant, and retrieve initial annotations like UniProt IDs.

#### Knowledge Graph Construction
The core of the system is a Neo4j graph database, which is populated by a series of modular Python scripts:

4.  **Base Graph Import (`1_neo4j_base_importer.py`):** Creates the foundational layer of the graph, consisting of `:Mutation` nodes from the GWAS data and the `:Gene` nodes they affect.
5.  **Ontology Integration (`2_go_importer.py`, `3_hpo_importer.py`):** Enriches the graph by importing the full Gene Ontology (GO) and Human Phenotype Ontology (HPO), creating `:GO_Term` and `:Phenotype` nodes and the rich hierarchical relationships between them (`IS_A`, `PART_OF`, etc.).
6.  **Pathway Integration (`4_reactome_importer.py`):** Adds `:Pathway` nodes from the Reactome database and connects them to the relevant `:Gene` nodes.
7.  **Clinical Significance (`5_clinvar_importer.py`):** Enriches the `:Mutation` nodes with clinical significance data (e.g., "Benign," "Pathogenic") and links them to `:Disease` nodes from the ClinVar database.

#### Unstructured Data Integration
8.  **Literature Retrieval (`6_pubmed_fetcher.py`):** Fetches relevant scientific abstracts from PubMed for the genes identified in the graph.
9.  **Knowledge Extraction (`7_ner_importer.py`):** Uses a GPU-accelerated `scispaCy` Named Entity Recognition (NER) model to read the downloaded abstracts, identify mentions of biomedical concepts (diseases, chemicals, etc.), and load them into the graph as generic `:Entity` nodes linked to the `:Paper` they were mentioned in.
10. **Entity Reconciliation (`8_ner_reconciliation.py`):** This crucial script links the unstructured and structured worlds. It intelligently merges the generic `:Entity` nodes created from the literature with the high-quality, "official" `:Gene`, `:Phenotype`, and `:Disease` nodes already in the graph.

#### RAG System Preparation
11. **Vector Embedding (`9_create_embeddings.py`):** Uses a PubMedBERT-based SentenceTransformer model to convert the semantic meaning of each PubMed abstract into a vector embedding. These embeddings are stored in a local ChromaDB vector database.

---

## How to Run the System

Follow these steps to build and query the knowledge graph from scratch.

### Step 0: Initial Setup

1.  **Environment:** Make sure your `gwas-env` conda environment is active.
    ```bash
    conda activate gwas-env
    ```
2.  **Dependencies:** Ensure all Python packages from our debugging sessions are installed.
3.  **Clean the Directory (Optional):** To start completely fresh, run the `clean.sh` script.
    ```bash
    bash clean.sh
    ```

### Step 1: Start the Database Service

Start the Neo4j Docker container in the background.
```bash
docker-compose up -d neo4j

Step 2: Run the Full Data Pipeline
Execute the master script. This will run all the data processing and import steps in the correct order.

bash run_pipeline.sh

Step 3: Query Your Knowledge Graph
Once the pipeline is complete, you can start the interactive RAG engine to ask questions.

python query_engine.py
