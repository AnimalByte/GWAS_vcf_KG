import os
import logging
from neo4j import GraphDatabase
import time

# --- Configuration ---
# Sets up logging to monitor the script's progress.
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Neo4j connection details
NEO4J_URI = os.getenv("NEO4J_URI", "bolt://localhost:7687")
NEO4J_USER = os.getenv("NEO4J_USER", "neo4j")
NEO4J_PASSWORD = os.getenv("NEO4J_PASSWORD", "password")


# --- Main Reconciliation Function ---
def reconcile_entities(driver):
    """
    Finds generic :Entity nodes created by the NER process and attempts to
    reconcile them with specific, typed nodes (:Gene, :Phenotype, :Disease)
    already in the graph by performing the matching logic in Python.
    """
    logging.info("--- Starting Entity Reconciliation Pipeline ---")

    with driver.session() as session:
        # --- Step 1: Fetch all necessary data from the graph ---
        
        logging.info("Fetching generic entities from the graph...")
        entity_result = session.run("MATCH (e:Entity) WHERE size(labels(e)) = 1 RETURN e.name AS name")
        generic_entities = {record['name'] for record in entity_result}
        logging.info(f"Found {len(generic_entities)} generic :Entity nodes to process.")
        
        logging.info("Fetching specific :Phenotype nodes from the graph...")
        phenotype_result = session.run("MATCH (p:Phenotype) RETURN p.name AS name")
        phenotype_map = {name.lower(): name for record in phenotype_result for name in (record['name'] if isinstance(record['name'], list) else [record['name']])}
        logging.info(f"Found {len(phenotype_map)} unique Phenotype names.")

        logging.info("Fetching specific :Disease nodes from the graph...")
        disease_result = session.run("MATCH (d:Disease) RETURN d.name AS name")
        disease_map = {name.lower(): name for record in disease_result for name in (record['name'] if isinstance(record['name'], list) else [record['name']])}
        logging.info(f"Found {len(disease_map)} unique Disease names.")
        
        logging.info("Fetching specific :Gene nodes from the graph...")
        gene_result = session.run("MATCH (g:Gene) RETURN g.symbol AS name")
        gene_map = {record['name'].lower(): record['name'] for record in gene_result}
        logging.info(f"Found {len(gene_map)} unique Gene names.")

        # --- Step 2: Perform reconciliation logic in Python ---

        logging.info("Finding matches between generic entities and specific nodes...")
        # Create separate lists for each type of reconciliation
        phenotype_matches = []
        disease_matches = []
        gene_matches = []

        for entity_name in generic_entities:
            lower_entity_name = entity_name.lower()
            if lower_entity_name in phenotype_map:
                phenotype_matches.append({"entity_name": entity_name, "specific_name": phenotype_map[lower_entity_name]})
            elif lower_entity_name in disease_map:
                disease_matches.append({"entity_name": entity_name, "specific_name": disease_map[lower_entity_name]})
            elif lower_entity_name in gene_map:
                gene_matches.append({"entity_name": entity_name, "specific_name": gene_map[lower_entity_name]})

        logging.info(f"Found {len(phenotype_matches)} Phenotype, {len(disease_matches)} Disease, and {len(gene_matches)} Gene reconciliations to perform.")

        # --- Step 3: Execute merge queries in Neo4j ---
        
        # **FIXED**: Using simpler, separate queries for each entity type.
        phenotype_merge_query = """
        UNWIND $rows AS row
        MATCH (generic:Entity {name: row.entity_name})
        MATCH (specific:Phenotype {name: row.specific_name})
        WITH generic, specific
        CALL apoc.refactor.mergeNodes([generic, specific], {properties: 'combine'}) YIELD node
        RETURN count(node) as merged_count
        """
        
        disease_merge_query = """
        UNWIND $rows AS row
        MATCH (generic:Entity {name: row.entity_name})
        MATCH (specific:Disease {name: row.specific_name})
        WITH generic, specific
        CALL apoc.refactor.mergeNodes([generic, specific], {properties: 'combine'}) YIELD node
        RETURN count(node) as merged_count
        """

        gene_merge_query = """
        UNWIND $rows AS row
        MATCH (generic:Entity {name: row.entity_name})
        MATCH (specific:Gene {symbol: row.specific_name})
        WITH generic, specific
        CALL apoc.refactor.mergeNodes([generic, specific], {properties: 'combine'}) YIELD node
        RETURN count(node) as merged_count
        """

        batch_size = 5000
        
        def run_merge(query, matches, label):
            total_merged = 0
            logging.info(f"Executing merge operations for : {label}...")
            for i in range(0, len(matches), batch_size):
                batch = matches[i:i + batch_size]
                result = session.run(query, rows=batch).data()
                merged_in_batch = sum(record['merged_count'] for record in result)
                total_merged += merged_in_batch
            logging.info(f"Total : {label} nodes merged: {total_merged}")

        run_merge(phenotype_merge_query, phenotype_matches, "Phenotype")
        run_merge(disease_merge_query, disease_matches, "Disease")
        run_merge(gene_merge_query, gene_matches, "Gene")


    logging.info(f"--- Entity Reconciliation Complete! ---")


# --- Main execution block ---
def main():
    """Orchestrates the NER pipeline."""
    try:
        driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASSWORD))
        
        reconcile_entities(driver)

        driver.close()

    except Exception as e:
        logging.error(f"An error occurred in the main pipeline: {e}")

if __name__ == "__main__":
    main()
