import os
import logging
import urllib.request
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

ENCODE_CCRE_URL = "https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.bed"

# --- Data Download Function ---
def download_file(url, directory, filename):
    """Downloads a file from a URL if it doesn't already exist."""
    if not os.path.exists(directory):
        os.makedirs(directory)
    filepath = os.path.join(directory, filename)
    if not os.path.exists(filepath):
        logging.info(f"Downloading {filename} from {url}...")
        urllib.request.urlretrieve(url, filepath)
        logging.info("Download complete.")
    else:
        logging.info(f"{filename} already exists. Skipping download.")
    return filepath

def import_encode_data(driver, bed_path):
    """
    Parses ENCODE cCRE data and loads it into the graph as :RegulatoryElement nodes.
    """
    if not os.path.exists(bed_path):
        logging.error(f"ENCODE data file not found at {bed_path}. Cannot proceed.")
        return False

    logging.info("Parsing ENCODE cCRE data...")
    try:
        df = pd.read_csv(bed_path, sep='\t', header=None, usecols=[0, 1, 2, 3, 4])
        df.columns = ['chrom', 'start', 'end', 'accession', 'label']
        
        records = df.to_dict('records')
        logging.info(f"Found {len(records)} candidate cis-Regulatory Elements (cCREs) to import.")

        query = """
        UNWIND $rows AS row
        CREATE (re:RegulatoryElement {
            accession: row.accession,
            label: row.label,
            location: row.chrom + ':' + row.start + '-' + row.end,
            chrom: row.chrom,
            start: toInteger(row.start),
            end: toInteger(row.end)
        })
        """

        with driver.session() as session:
            session.run("CREATE INDEX IF NOT EXISTS FOR (n:RegulatoryElement) ON (n.accession)").consume()
            session.run("CREATE INDEX IF NOT EXISTS FOR (n:RegulatoryElement) ON (n.chrom)").consume()
            session.run("CREATE RANGE INDEX IF NOT EXISTS FOR (n:RegulatoryElement) ON (n.start, n.end)").consume()
            
            logging.info("Importing ENCODE cCREs into Neo4j...")
            batch_size = 50000
            for i in range(0, len(records), batch_size):
                batch = records[i:i + batch_size]
                summary = session.run(query, rows=batch).consume()
                logging.info(f"Imported chunk {i//batch_size + 1}. Nodes created: {summary.counters.nodes_created}")

        logging.info("ENCODE data import complete.")
        return True

    except Exception as e:
        logging.error(f"An error occurred during ENCODE import: {e}")
        return False

def link_variants_to_regulatory_elements(driver):
    """
    Creates :IN_REGULATORY_ELEMENT relationships between :Mutation nodes
    and the :RegulatoryElement nodes they overlap with. This version is optimized
    and provides progress updates.
    """
    logging.info("Linking GWAS variants to overlapping ENCODE regulatory elements...")
    try:
        with driver.session() as session:
            # 1. Fetch all mutations from the graph into memory
            logging.info("Fetching all mutation data from graph...")
            mutations_result = session.run("MATCH (m:Mutation) WHERE m.chrom IS NOT NULL AND m.pos IS NOT NULL RETURN m.id AS id, m.chrom AS chrom, m.pos AS pos")
            mutations = [{"id": record["id"], "chrom": record["chrom"], "pos": int(record["pos"])} for record in mutations_result]
            total_mutations = len(mutations)
            
            if total_mutations == 0:
                logging.warning("No mutations found in the graph to link. Skipping.")
                return

            logging.info(f"Found {total_mutations} mutations to process.")
            
            # --- Extensive Debugging ---
            print("\n--- DEBUGGING: Chromosome Name Mismatch Check ---")
            sample_mutation_chroms = {m['chrom'] for m in mutations[:5]}
            print(f"Sample of chromosome names from Mutation nodes: {sample_mutation_chroms}")
            
            re_chroms_result = session.run("MATCH (re:RegulatoryElement) RETURN DISTINCT re.chrom LIMIT 5").value()
            print(f"Sample of chromosome names from RegulatoryElement nodes: {re_chroms_result}")
            print("--------------------------------------------------\n")
            # --- End Debugging ---

            # 2. Process in batches to provide progress feedback
            # *** CORRECTION: Handle both 'chr1' and '1' style chromosome names ***
            link_query = """
            UNWIND $batch AS mutation
            MATCH (m:Mutation {id: mutation.id})
            MATCH (re:RegulatoryElement)
            WHERE (re.chrom = mutation.chrom OR re.chrom = 'chr' + mutation.chrom)
              AND mutation.pos >= re.start 
              AND mutation.pos <= re.end
            MERGE (m)-[:IN_REGULATORY_ELEMENT]->(re)
            """
            
            batch_size = 5000  # Smaller batch size for this more complex query
            total_rels_created = 0
            for i in range(0, total_mutations, batch_size):
                batch = mutations[i:i + batch_size]
                summary = session.run(link_query, batch=batch).consume()
                rels_created = summary.counters.relationships_created
                total_rels_created += rels_created
                logging.info(f"Processed batch {i // batch_size + 1}/{ -(-total_mutations // batch_size) }... {i + len(batch)}/{total_mutations} mutations processed. New relationships: {rels_created}")

            logging.info(f"Finished linking variants. Total relationships created: {total_rels_created}")
            
    except Exception as e:
        logging.error(f"An error occurred while linking variants to regulatory elements: {e}")

# --- Main execution block ---
def main_encode_pipeline():
    """Orchestrates the download, import, and linking of ENCODE cCRE data."""
    try:
        logging.info("--- Starting ENCODE cCRE Data Import Pipeline ---")
        time.sleep(5)
        
        download_dir = "downloads"
        data_path = download_file(ENCODE_CCRE_URL, download_dir, "GRCh38-cCREs.bed")

        driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASSWORD))
        
        # Delete existing nodes in batches to avoid out-of-memory errors.
        with driver.session() as session:
            logging.info("Clearing any previous ENCODE data...")
            while True:
                result = session.run("MATCH (re:RegulatoryElement) WITH re LIMIT 100000 DETACH DELETE re RETURN count(re)")
                count = result.single()[0]
                if count == 0:
                    break
                logging.info(f"Deleted {count} old ENCODE nodes.")
            logging.info("Finished clearing old ENCODE data.")
            
        if import_encode_data(driver, data_path):
            link_variants_to_regulatory_elements(driver)
        
        driver.close()
        logging.info("--- ENCODE Data Import Pipeline Finished Successfully! ---")

    except Exception as e:
        logging.error(f"An error occurred in the ENCODE pipeline: {e}")

if __name__ == "__main__":
    main_encode_pipeline()
