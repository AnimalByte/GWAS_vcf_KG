import os
import logging
import urllib.request
import gzip
import shutil
import pandas as pd
from neo4j import GraphDatabase
import time

# --- Configuration ---
# Sets up logging to monitor the script's progress.
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Neo4j connection details
NEO4J_URI = os.getenv("NEO4J_URI", "bolt://localhost:7687")
NEO4J_USER = os.getenv("NEO4J_USER", "neo4j")
NEO4J_PASSWORD = os.getenv("NEO4J_PASSWORD", "password")

# STRING-DB Configuration
STRING_VERSION = "12.0" 
CONFIDENCE_THRESHOLD = 700  # Medium-high confidence (scale is 1-1000)

# --- Data Download Function ---
def download_and_decompress(url, directory, filename):
    """Downloads and decompresses a .gz file if the uncompressed version doesn't already exist."""
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    gz_path = os.path.join(directory, filename)
    uncompressed_path = gz_path.replace('.gz', '')

    if not os.path.exists(uncompressed_path):
        if not os.path.exists(gz_path):
            logging.info(f"Downloading {filename} from {url}...")
            urllib.request.urlretrieve(url, gz_path)
            logging.info("Download complete.")
        
        logging.info(f"Decompressing {filename}...")
        with gzip.open(gz_path, 'rb') as f_in:
            with open(uncompressed_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        logging.info(f"Decompression complete: {uncompressed_path}")
    else:
        logging.info(f"{uncompressed_path} already exists. Skipping download and decompression.")
        
    return uncompressed_path

def import_ppi_data(driver, links_path, alias_path):
    """
    Parses STRING-DB data to create protein-protein interaction relationships
    between existing :Gene nodes in the graph. It will create new :Gene nodes if they don't exist.
    """
    if not all([os.path.exists(links_path), os.path.exists(alias_path)]):
        logging.error("One or more STRING data files not found. Cannot proceed.")
        return

    # 1. Create a mapping from STRING ID to Gene Symbol using a more robust method
    logging.info("Creating map from STRING Protein ID to Gene Symbol...")
    try:
        # Assign column names manually as the file has no header
        column_names = ['string_protein_id', 'alias', 'source']
        # Use comment='#' to correctly skip the header line
        df_aliases = pd.read_csv(alias_path, sep='\t', header=None, names=column_names, comment='#', low_memory=False)
        
        # Prioritize 'Ensembl_HGNC_symbol' source for the highest quality mapping
        hgnc_aliases = df_aliases[df_aliases['source'].str.contains("Ensembl_HGNC_symbol", na=False)]
        
        string_to_symbol = pd.Series(hgnc_aliases.alias.values, index=hgnc_aliases.string_protein_id).to_dict()
        logging.info(f"Mapped {len(string_to_symbol)} STRING IDs using primary HGNC symbols.")
        
        # For STRING IDs not yet mapped, use the first available alias as a fallback
        remaining_aliases = df_aliases[~df_aliases['string_protein_id'].isin(string_to_symbol.keys())]
        fallback_map = remaining_aliases.groupby('string_protein_id')['alias'].first().to_dict()
        
        string_to_symbol.update(fallback_map)
        logging.info(f"Mapped a total of {len(string_to_symbol)} STRING IDs after fallback.")

    except Exception as e:
        logging.error(f"Failed to parse alias file: {e}")
        return

    # 2. Process interaction links
    logging.info(f"Processing interaction data and filtering for confidence > {CONFIDENCE_THRESHOLD}...")
    interactions = []
    with open(links_path, 'r', encoding='utf-8') as f:
        next(f) # Skip header
        for line in f:
            parts = line.strip().split()
            protein1, protein2, score = parts[0], parts[1], int(parts[-1]) # Combined score is the last column
            
            if score >= CONFIDENCE_THRESHOLD:
                gene1 = string_to_symbol.get(protein1)
                gene2 = string_to_symbol.get(protein2)
                
                if gene1 and gene2 and gene1 != gene2:
                    if gene1 > gene2:
                        gene1, gene2 = gene2, gene1
                    
                    interactions.append({
                        "gene1": gene1,
                        "gene2": gene2,
                        "score": score
                    })

    # Remove duplicate interactions
    if interactions:
        interactions_df = pd.DataFrame(interactions).drop_duplicates()
        interactions = interactions_df.to_dict('records')

    logging.info(f"Found {len(interactions)} unique, high-confidence interactions to import.")
    
    if not interactions:
        logging.warning("No interactions to import after filtering.")
        return

    # 3. Import into Neo4j
    # *** CORRECTION: Use MERGE for genes to create them if they don't exist ***
    query = """
    UNWIND $rows AS row
    MERGE (g1:Gene {symbol: row.gene1})
    MERGE (g2:Gene {symbol: row.gene2})
    MERGE (g1)-[r:INTERACTS_WITH {source: 'STRING-DB'}]-(g2)
    ON CREATE SET r.score = row.score
    ON MATCH SET r.score = CASE WHEN r.score < row.score THEN row.score ELSE r.score END
    """

    with driver.session() as session:
        logging.info("Importing PPI relationships into Neo4j...")
        batch_size = 50000
        total_rels_created = 0
        total_nodes_created = 0

        for i in range(0, len(interactions), batch_size):
            batch = interactions[i:i + batch_size]
            result = session.run(query, rows=batch)
            summary = result.consume()
            
            rels = summary.counters.relationships_created
            nodes = summary.counters.nodes_created
            total_rels_created += rels
            total_nodes_created += nodes

            logging.info(f"Imported chunk {i // batch_size + 1}. Nodes created: {nodes}, Relationships created: {rels}")

        logging.info(f"PPI data import complete. Total nodes created: {total_nodes_created}, Total relationships created: {total_rels_created}.")

# --- Main execution block ---
def main_ppi_pipeline():
    """Orchestrates the download and import of STRING-DB PPI data."""
    try:
        logging.info("--- Starting STRING-DB Protein-Protein Interaction Import Pipeline ---")
        time.sleep(5)
        
        download_dir = "downloads"
        
        links_url = f"https://stringdb-static.org/download/protein.links.full.v{STRING_VERSION}/9606.protein.links.full.v{STRING_VERSION}.txt.gz"
        alias_url = f"https://stringdb-static.org/download/protein.aliases.v{STRING_VERSION}/9606.protein.aliases.v{STRING_VERSION}.txt.gz"
        
        links_file_path = download_and_decompress(links_url, download_dir, f"9606.protein.links.full.v{STRING_VERSION}.txt.gz")
        alias_file_path = download_and_decompress(alias_url, download_dir, f"9606.protein.aliases.v{STRING_VERSION}.txt.gz")

        driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASSWORD))
        
        # Clear any previous STRING-DB PPI data
        with driver.session() as session:
            logging.info("Clearing any previous STRING-DB PPI relationships...")
            session.run("MATCH ()-[r:INTERACTS_WITH {source: 'STRING-DB'}]-() DELETE r").consume()

        import_ppi_data(driver, links_file_path, alias_file_path)
        
        driver.close()
        logging.info("--- STRING-DB Import Pipeline Finished Successfully! ---")

    except Exception as e:
        logging.error(f"An error occurred in the PPI pipeline: {e}")

if __name__ == "__main__":
    main_ppi_pipeline()
