import os
import time
import logging
from neo4j import GraphDatabase
from Bio import Entrez, Medline
from datetime import datetime, timedelta

# --- Configuration ---
# Sets up logging to monitor the script's progress.
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Neo4j connection details
NEO4J_URI = os.getenv("NEO4J_URI", "bolt://localhost:7687")
NEO4J_USER = os.getenv("NEO4J_USER", "neo4j")
NEO4J_PASSWORD = os.getenv("NEO4J_PASSWORD", "password")

# NCBI Entrez API configuration
Entrez.email = "your_email@example.com"  # Please change this to your email
MAX_PAPERS_PER_GENE = 15  # Max number of abstracts to download per gene
SEARCH_YEARS_PRIMARY = 7 # *** The primary, most recent search window ***
SEARCH_YEARS_FALLBACK = 15 # *** The fallback window if no recent results are found ***

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

def search_pubmed(gene_symbol, years):
    """Performs a PubMed search for a given gene and time window."""
    end_date = datetime.now()
    start_date = end_date - timedelta(days=years * 365)
    start_date_str = start_date.strftime("%Y/%m/%d")
    end_date_str = end_date.strftime("%Y/%m/%d")
    
    search_term = (
        f'({gene_symbol}[Gene Name] OR "{gene_symbol}"[All Fields]) AND '
        f'"human"[Organism] AND hasabstract[Filter] AND '
        f'("{start_date_str}"[Date - Publication] : "{end_date_str}"[Date - Publication])'
    )
    
    logging.info(f"Searching PubMed for '{gene_symbol}' within the last {years} years...")
    handle = Entrez.esearch(db="pubmed", term=search_term, retmax=MAX_PAPERS_PER_GENE)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]

def fetch_pubmed_abstracts(gene_symbol):
    """
    Searches PubMed for a given gene symbol and downloads abstracts.
    It first tries a 7-year window and falls back to a 15-year window if needed.
    """
    try:
        # --- Primary search attempt (last 7 years) ---
        pmids = search_pubmed(gene_symbol, SEARCH_YEARS_PRIMARY)

        # --- Fallback search attempt (last 15 years) ---
        if not pmids:
            logging.warning(f"No results found for '{gene_symbol}' in the last {SEARCH_YEARS_PRIMARY} years. Expanding search to {SEARCH_YEARS_FALLBACK} years.")
            pmids = search_pubmed(gene_symbol, SEARCH_YEARS_FALLBACK)

        if not pmids:
            logging.info(f"No results found for '{gene_symbol}' even in the last {SEARCH_YEARS_FALLBACK} years. Skipping.")
            return

        logging.info(f"Found {len(pmids)} PMIDs for '{gene_symbol}'.")

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
