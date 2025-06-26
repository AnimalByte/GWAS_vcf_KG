import os
import time
import logging
from neo4j import GraphDatabase
from Bio import Entrez, Medline

# --- Configuration ---
# Sets up logging to monitor the script's progress.
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Neo4j connection details
NEO4J_URI = os.getenv("NEO4J_URI", "bolt://localhost:7687")
NEO4J_USER = os.getenv("NEO4J_USER", "neo4j")
NEO4J_PASSWORD = os.getenv("NEO4J_PASSWORD", "password")

# NCBI Entrez API configuration
# It's polite to tell NCBI who you are
Entrez.email = "dannysoriano89@gmail.com"  # Please change this to your email
MAX_PAPERS_PER_GENE = 15  # Max number of abstracts to download per gene

# Directory to save the downloaded abstracts
OUTPUT_DIR = "pubmed_abstracts"

def get_genes_from_neo4j(driver):
    """Queries the Neo4j database to get a list of all gene symbols."""
    logging.info("Fetching gene symbols from Neo4j...")
    genes = []
    try:
        with driver.session() as session:
            # This query fetches all unique gene symbols from your graph
            result = session.run("MATCH (g:Gene) WHERE g.symbol IS NOT NULL RETURN DISTINCT g.symbol AS symbol")
            genes = [record["symbol"] for record in result]
        logging.info(f"Found {len(genes)} unique gene symbols in the graph to process.")
        return genes
    except Exception as e:
        logging.error(f"Failed to query genes from Neo4j: {e}")
        return []

def fetch_pubmed_abstracts(gene_symbol):
    """
    Searches PubMed for a given gene symbol and downloads abstracts.
    """
    # This search term is robust, looking for the official gene name or the symbol in any field.
    search_term = f'({gene_symbol}[Gene Name] OR "{gene_symbol}"[All Fields]) AND "human"[Organism] AND hasabstract[Filter]'
    logging.info(f"Searching PubMed for gene '{gene_symbol}'...")
    
    try:
        # Search PubMed for article IDs (PMIDs)
        handle = Entrez.esearch(db="pubmed", term=search_term, retmax=MAX_PAPERS_PER_GENE)
        record = Entrez.read(handle)
        handle.close()
        pmids = record["IdList"]

        logging.info(f"Found {len(pmids)} PMIDs for '{gene_symbol}'.")

        if not pmids:
            return

        # Fetch the actual records for the found PMIDs
        handle = Entrez.efetch(db="pubmed", id=pmids, rettype="medline", retmode="text")
        records = Medline.parse(handle)
        
        # Create a directory for the gene
        gene_dir = os.path.join(OUTPUT_DIR, gene_symbol)
        if not os.path.exists(gene_dir):
            os.makedirs(gene_dir)

        # Save each abstract to a file
        saved_count = 0
        for rec in records:
            abstract = rec.get("AB")
            pmid = rec.get("PMID")
            if abstract:
                filepath = os.path.join(gene_dir, f"{pmid}.txt")
                with open(filepath, "w", encoding="utf-8") as f:
                    f.write(abstract)
                saved_count += 1
        
        logging.info(f"Successfully saved {saved_count} abstracts for {gene_symbol}.")

        # Be polite to the NCBI API and don't spam requests
        time.sleep(1)

    except Exception as e:
        logging.error(f"An error occurred while fetching data for {gene_symbol}: {e}")


def main():
    """Main function to run the PubMed data fetching pipeline."""
    try:
        driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASSWORD))
        
        genes = get_genes_from_neo4j(driver)
        
        if not genes:
            logging.warning("No genes found in Neo4j. Cannot fetch PubMed data.")
            return

        # Create the main output directory if it doesn't exist
        if not os.path.exists(OUTPUT_DIR):
            os.makedirs(OUTPUT_DIR)
            
        logging.info(f"Starting PubMed abstract download for {len(genes)} genes...")
        for gene in genes:
            fetch_pubmed_abstracts(gene)
            
        driver.close()
        logging.info("PubMed abstract fetching complete!")

    except Exception as e:
        logging.error(f"An error occurred in the main pipeline: {e}")

if __name__ == "__main__":
    main()
