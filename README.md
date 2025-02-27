# README.md

# Annotating VCF Files with GO Ontologies Using VEP and Neo4j

This project processes a public dataset, extracts significant genetic variants, annotates them with GO ontologies using VEP, and integrates the results into a Neo4j knowledge graph.

---

## 📂 Files Overview

### `demonstration.sh`
This Bash script automates the process of extracting genetic variant information from a VCF file, filtering significant variants, converting them back to VCF format, and annotating them with GO ontology terms using VEP.

#### 🔹 Steps:
1. **Extract GWAS Summary Statistics**  
   - Uses `bcftools` to extract information from a VCF file.  
   - Outputs `data/gwas_summary_stats.txt`.

2. **Filter Significant Variants**  
   - Filters variants where `LP > -log10(5e-8)`.  
   - Outputs `results/significant_variants.txt`.

3. **Convert Significant Variants to VCF Format**  
   - Reformats filtered variants back to a VCF file.  
   - Outputs `results/significant_variants.vcf`.

4. **Annotate Significant Variants with VEP**  
   - Uses Ensembl's `VEP` to annotate variants with gene ontology (GO) terms.  
   - Outputs `results/vep_annotated.txt`.

5. **Map Variants to Genes**  
   - Extracts variant-gene relationships from VEP results.  
   - Outputs `annotations/snp_gene_map.txt`.

---

### `pathway_extraction.py`
This Python script processes the annotated variant data and integrates it into a Neo4j knowledge graph.

#### 🔹 Steps:
1. **Import GO Ontology from OWL File**  
   - Loads GO terms into Neo4j using Neosemantics (`n10s`).

2. **Parse and Import Gene-GO Annotations from GAF File**  
   - Parses a GAF annotation file.  
   - Links genes to GO terms in Neo4j.

3. **Import VEP Annotated Variants into Neo4j**  
   - Reads the VEP-annotated variant file.  
   - Links genetic variants to their associated genes in Neo4j.

#### 🔧 Requirements:
- Python packages: `pandas`, `neo4j`, `logging`
- A running Neo4j instance with the Neosemantics plugin (`n10s`)

---

## 🚀 Running the Pipeline

1. **Run the Bash script to process the VCF file:**
   ```sh
   chmod +x demonstration.sh
   ./demonstration.sh
2. ** Run the Python script to import the data into Neo4j:
   python pathway_extraction.py

