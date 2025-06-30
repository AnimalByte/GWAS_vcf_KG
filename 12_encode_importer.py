import os
import logging
import urllib.request
import pandas as pd
from neo4j import GraphDatabase
import time
import gzip
import shutil

# --- Configuration ---
# Sets up logging to monitor the script's progress.
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Neo4j connection details
NEO4J_URI = os.getenv("NEO4J_URI", "bolt://localhost:7687")
NEO4J_USER = os.getenv("NEO4J_USER", "neo4j")
NEO4J_PASSWORD = os.getenv("NEO4J_PASSWORD", "password")

# --- Data URLs ---
ENCODE_CCRE_URL = "https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.bed"
# CORRECTED URL for the cell type activity data matrix
ENCODE_CELL_TYPE_URL = "https://media.wenglab.org/screen/human/V4/GRCh38-cCREs/GRCh38-cCREs.DNase-Z-score.tsv.gz"
ACTIVITY_THRESHOLD = 1.64  # Standard Z-score for significant activity

# --- Data Download Function ---
def download_file(url, directory, filename, force_decompress=False):
    """Downloads a file and optionally decompresses it."""
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    filepath = os.path.join(directory, filename)
    
    if not os.path.exists(filepath):
        logging.info(f"Downloading {filename} from {url}...")
        urllib.request.urlretrieve(url, filepath)
        logging.info("Download complete.")
    else:
        logging.info(f"{filename} already exists. Skipping download.")

    if filename.endswith(".gz"):
        uncompressed_path = filepath[:-3]
        if not os.path.exists(uncompressed_path) or force_decompress:
            logging.info(f"Decompressing {filename}...")
            with gzip.open(filepath, 'rb') as f_in:
                with open(uncompressed_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            logging.info(f"Decompression complete: {uncompressed_path}")
        return uncompressed_path
    
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

def add_cell_type_specificity(driver, cell_type_path):
    """
    Reads the cell-type activity matrix and updates the :RegulatoryElement nodes
    with a list of biosamples where they are active.
    """
    if not os.path.exists(cell_type_path):
        logging.error(f"Cell type activity file not found at {cell_type_path}. Skipping this step.")
        return

    logging.info(f"Processing cell type activity data from {cell_type_path}...")
    
    update_query = """
    UNWIND $rows as row
    MATCH (re:RegulatoryElement {accession: row.accession})
    SET re.active_in = row.active_cell_types
    """

    try:
        # Process the large file in chunks to manage memory
        chunk_size = 25000 
        reader = pd.read_csv(cell_type_path, sep='\t', index_col=0, chunksize=chunk_size)
        total_processed = 0

        with driver.session() as session:
            for i, chunk in enumerate(reader):
                records_to_update = []
                # For each row (cCRE) in the chunk...
                for accession, data in chunk.iterrows():
                    # Find all columns (cell types) where the Z-score is above our threshold
                    active_cells = data[data > ACTIVITY_THRESHOLD].index.tolist()
                    if active_cells:
                        records_to_update.append({
                            "accession": accession,
                            "active_cell_types": active_cells
                        })
                
                # If we found any active elements, update them in the database
                if records_to_update:
                    session.run(update_query, rows=records_to_update).consume()
                
                total_processed += len(chunk)
                logging.info(f"Processed chunk {i+1}. Total cCREs processed: {total_processed}")

        logging.info("Successfully added cell type specificity to regulatory element nodes.")
    except Exception as e:
        logging.error(f"Failed to process cell type specificity data: {e}", exc_info=True)


def link_variants_to_regulatory_elements(driver):
    """
    Creates :IN_REGULATORY_ELEMENT relationships between :Mutation nodes
    and the :RegulatoryElement nodes they overlap with.
    """
    logging.info("Linking GWAS variants to overlapping ENCODE regulatory elements...")
    try:
        with driver.session() as session:
            link_query = """
            CALL apoc.periodic.iterate(
                "MATCH (m:Mutation) WHERE m.chrom IS NOT NULL AND m.pos IS NOT NULL RETURN m",
                "WITH m
                 MATCH (re:RegulatoryElement)
                 WHERE (re.chrom = m.chrom OR re.chrom = 'chr' + m.chrom OR m.chrom = 'chr' + re.chrom)
                   AND toInteger(m.pos) >= re.start
                   AND toInteger(m.pos) <= re.end
                 MERGE (m)-[:IN_REGULATORY_ELEMENT]->(re)",
                {batchSize: 5000, parallel: true, retries: 3}
            ) YIELD batches, total, timeTaken
            RETURN batches, total, timeTaken
            """
            
            logging.info("Executing high-performance linking query with apoc.periodic.iterate...")
            result = session.run(link_query).single()
            if result:
                logging.info(f"Finished linking variants. Processed {result['total']} mutations in {result['batches']} batches. Time taken: {result['timeTaken']}s.")
            else:
                logging.error("Linking query did not return expected results. Check APOC installation and graph data.")

    except Exception as e:
        logging.error(f"An error occurred while linking variants: {e}", exc_info=True)
        logging.warning("APOC-based linking failed. You may need to install the APOC plugin.")


# --- Main execution block ---
def main_encode_pipeline():
    """Orchestrates the download, import, and linking of ENCODE cCRE data."""
    try:
        logging.info("--- Starting ENCODE cCRE Data Import Pipeline ---")
        time.sleep(5)
        
        download_dir = "downloads"
        ccre_bed_path = download_file(ENCODE_CCRE_URL, download_dir, "GRCh38-cCREs.bed")
        cell_type_path = download_file(ENCODE_CELL_TYPE_URL, download_dir, "GRCh38-cCREs.DNase-Z-score.tsv.gz")

        driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASSWORD))
        
        with driver.session() as session:
            logging.info("Clearing any previous ENCODE data...")
            while True:
                result = session.run("MATCH (re:RegulatoryElement) WITH re LIMIT 100000 DETACH DELETE re RETURN count(re)").consume()
                if result.counters.nodes_deleted == 0:
                    break
                logging.info(f"Deleted a batch of old ENCODE nodes.")
            logging.info("Finished clearing old ENCODE data.")
            
        if import_encode_data(driver, ccre_bed_path):
            add_cell_type_specificity(driver, cell_type_path)
            link_variants_to_regulatory_elements(driver)
        
        driver.close()
        logging.info("--- ENCODE Data Import Pipeline Finished Successfully! ---")

    except Exception as e:
        logging.error(f"An error occurred in the ENCODE pipeline: {e}", exc_info=True)

if __name__ == "__main__":
    main_encode_pipeline()
