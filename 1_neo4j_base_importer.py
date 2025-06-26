import os
import logging
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

# --- Path to the new, correctly annotated data file ---
ANNOTATED_TSV_PATH = "results/egIwc7NRt4hou5yo.txt"

def import_base_graph_data(driver, tsv_path):
    """
    Parses the annotated TSV file and loads the base graph of Mutations,
    Genes, and their relationships into Neo4j. This version reads the file
    line-by-line with extensive debugging to be robust against formatting errors.
    """
    if not os.path.exists(tsv_path):
        logging.error(f"Annotated TSV file not found at {tsv_path}. Cannot proceed.")
        return

    logging.info(f"Step 1: Parsing annotated TSV file line-by-line: {tsv_path}")

    header = []
    all_records = []
    skipped_lines = 0

    try:
        with open(tsv_path, 'r', encoding='utf-8') as f:
            for i, line in enumerate(f):
                # Find and process the header line
                if line.startswith('#Uploaded_variation'):
                    header = [h.strip() for h in line.strip().split('\t')]
                    # Clean the first column name
                    if header and header[0] == '#Uploaded_variation':
                        header[0] = 'Uploaded_variation'
                    logging.info(f"DEBUG: Found header with {len(header)} columns: {header}")
                    continue
                # Skip other comment lines or empty lines
                if line.startswith('#') or not line.strip():
                    continue

                # Split the data line by tabs
                fields = [field.strip() for field in line.strip().split('\t')]

                # **DEBUGGING**: Check for column count consistency to handle malformed lines
                if len(fields) == len(header):
                    try:
                        record_dict = dict(zip(header, fields))
                        all_records.append(record_dict)
                    except Exception as e:
                        logging.error(f"Error creating dictionary on line {i+1}: {e}")
                        skipped_lines += 1
                else:
                    logging.warning(f"Skipping malformed line {i+1}. Expected {len(header)} columns, found {len(fields)}.")
                    # **DEBUGGING**: Print the content of the problematic line
                    logging.warning(f"  - Content: {line.strip()}")
                    skipped_lines += 1

        if not header:
            logging.error("Could not find a valid header row in the file.")
            return

        logging.info(f"File parsing summary: Total lines read: {i+1}, Successfully parsed: {len(all_records)}, Skipped lines: {skipped_lines}")

        # --- DEBUGGING: Show a sample of the successfully parsed data ---
        if all_records:
            logging.info("Sample of successfully parsed records:")
            # Use pandas just for pretty printing the sample
            print(pd.DataFrame(all_records[:3]).to_string())
        # --- END DEBUGGING ---


        # This Cypher query creates the core of the knowledge graph
        query = """
        UNWIND $rows AS row
        // Create the Mutation node using its unique ID from the VEP output
        MERGE (m:Mutation {id: row.Uploaded_variation})
        SET m.location = row.Location,
            m.chrom = split(row.Location, ':')[0],
            m.pos = split(split(row.Location, ':')[1], '-')[0],
            m.ref = row.REF_ALLELE,
            m.alt = row.Allele,
            m.consequence = row.Consequence,
            m.impact = row.IMPACT

        // Create the Gene node if it doesn't exist and link the Mutation to it
        WITH m, row
        WHERE row.SYMBOL IS NOT NULL AND row.SYMBOL <> '-' AND row.SYMBOL <> ''
        MERGE (g:Gene {symbol: row.SYMBOL})
        MERGE (m)-[:AFFECTS]->(g)
        """

        with driver.session() as session:
            # Create indexes for faster merging and lookups
            session.run("CREATE INDEX mutation_id_index IF NOT EXISTS FOR (n:Mutation) ON (n.id)").consume()
            session.run("CREATE INDEX gene_symbol_index IF NOT EXISTS FOR (n:Gene) ON (n.symbol)").consume()

            logging.info("Step 2: Importing Mutations and Genes into Neo4j...")

            chunk_size = 50000
            total_nodes_created = 0
            total_rels_created = 0

            for i in range(0, len(all_records), chunk_size):
                batch = all_records[i:i+chunk_size]
                result = session.run(query, rows=batch)
                summary = result.consume()
                nodes_created = summary.counters.nodes_created
                rels_created = summary.counters.relationships_created
                total_nodes_created += nodes_created
                total_rels_created += rels_created
                logging.info(f"Imported chunk {i // chunk_size + 1}. Nodes created: {nodes_created}, Relationships created: {rels_created}")

        logging.info(f"Base graph import complete. Total nodes created: {total_nodes_created}, Total relationships created: {total_rels_created}.")

    except FileNotFoundError:
        logging.error(f"Annotated file not found at {tsv_path}. Please check the path.")
    except Exception as e:
        logging.error(f"An error occurred during the import process: {e}")

# --- Main execution block ---
def main_pipeline():
    """Orchestrates the base graph import."""
    try:
        logging.info("--- Starting Neo4j Base Import Pipeline ---")

        time.sleep(10) # Give Neo4j a moment to be ready

        driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASSWORD))

        import_base_graph_data(driver, ANNOTATED_TSV_PATH)

        driver.close()
        logging.info("--- Neo4j Base Import Pipeline Finished Successfully! ---")

    except Exception as e:
        logging.error(f"An error occurred in the main pipeline: {e}")

if __name__ == "__main__":
    main_pipeline()
