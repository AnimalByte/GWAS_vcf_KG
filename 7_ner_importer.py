import os
import logging
import spacy
from neo4j import GraphDatabase
import time

# --- Configuration ---
# Sets up logging to monitor the script's progress.
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Neo4j connection details
NEO4J_URI = os.getenv("NEO4J_URI", "bolt://localhost:7687")
NEO4J_USER = os.getenv("NEO4J_USER", "neo4j")
NEO4J_PASSWORD = os.getenv("NEO4J_PASSWORD", "password")

# --- Paths and Model Configuration ---
# Directory where the downloaded abstracts are stored.
ABSTRACTS_DIR = "pubmed_abstracts"
# Using the largest and most accurate scispaCy model.
SPACY_MODEL_NAME = "en_core_sci_lg"

# --- Main NER and Import Function ---
def extract_and_load_entities(driver):
    """
    Processes downloaded PubMed abstracts, extracts biomedical entities using a
    scispaCy model, and loads them into the Neo4j graph.
    """
    logging.info("--- Starting NER and Graph Enrichment Pipeline ---")

    # --- Step 1: Load the scispaCy NER Model ---
    logging.info("Checking for GPU and attempting to enable it for spaCy...")
    if spacy.prefer_gpu():
        logging.info("GPU enabled for spaCy!")
    else:
        logging.warning("GPU not available for spaCy, using CPU. This may be slower.")

    logging.info(f"Loading NER model: {SPACY_MODEL_NAME}...")
    try:
        nlp = spacy.load(SPACY_MODEL_NAME)
        logging.info("NER model loaded successfully.")
    except OSError:
        logging.error(f"Model '{SPACY_MODEL_NAME}' not found. Please run the correct installation commands.")
        logging.error("pip install scispacy==0.5.1")
        logging.error("pip install https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.1/en_core_sci_lg-0.5.1.tar.gz")
        return

    # --- Step 2: Prepare Neo4j Indexes ---
    with driver.session() as session:
        logging.info("Creating indexes for new entity types...")
        session.run("CREATE INDEX paper_pmid_index IF NOT EXISTS FOR (n:Paper) ON (n.pmid)").consume()
        session.run("CREATE INDEX entity_name_index IF NOT EXISTS FOR (n:Entity) ON (n.name)").consume()


    # --- Step 3: Process Abstracts and Load into Graph ---
    logging.info(f"Scanning for abstracts in: '{ABSTRACTS_DIR}'")
    
    # This query uses APOC to dynamically add labels to the nodes.
    query = """
    MERGE (p:Paper {pmid: $pmid})
    SET p.gene_context = $gene_symbol
    
    WITH p
    UNWIND $entities AS entity_data
    
    MERGE (e:Entity {name: entity_data.text})
    // This WITH is crucial for chaining the operations correctly.
    WITH p, e, entity_data
    CALL apoc.create.addLabels(e, [entity_data.label]) YIELD node
    
    MERGE (node)-[:MENTIONED_IN]->(p)
    """
    
    total_processed = 0
    # Walk through the directory structure.
    for gene_symbol in os.listdir(ABSTRACTS_DIR):
        gene_dir = os.path.join(ABSTRACTS_DIR, gene_symbol)
        if os.path.isdir(gene_dir):
            for filename in os.listdir(gene_dir):
                if filename.endswith(".txt"):
                    filepath = os.path.join(gene_dir, filename)
                    pmid = os.path.splitext(filename)[0]

                    with open(filepath, 'r', encoding='utf-8') as f:
                        abstract_text = f.read()
                    
                    doc = nlp(abstract_text)
                    
                    # We will accept any entity the model finds.
                    entities_to_load = []
                    for ent in doc.ents:
                         if len(ent.text) > 3: # Skip very short, likely non-specific entities
                            entities_to_load.append({
                                "text": ent.text.lower(), 
                                "label": ent.label_.capitalize()
                            })

                    if entities_to_load:
                        # Remove duplicate entities for the same paper before loading
                        unique_entities = [dict(t) for t in {tuple(d.items()) for d in entities_to_load}]
                        with driver.session() as session:
                            session.run(query, pmid=pmid, gene_symbol=gene_symbol, entities=unique_entities)
                    
                    total_processed += 1
                    if total_processed % 100 == 0:
                        logging.info(f"Processed {total_processed} abstracts...")

    logging.info(f"--- NER processing complete! Processed {total_processed} total abstracts. ---")

    # --- Step 4: Final Verification ---
    with driver.session() as session:
        result = session.run("MATCH (p:Paper) RETURN count(p) as count").single()
        paper_count = result['count'] if result else 0
        logging.info(f"VERIFICATION: Found {paper_count} :Paper nodes in the database.")
        if paper_count > 0:
            logging.info("SUCCESS: Data was loaded into the graph.")
        else:
            logging.error("FAILURE: No :Paper nodes were created.")


# --- Main execution block ---
def main():
    """Orchestrates the NER pipeline."""
    try:
        driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASSWORD))
        
        # Clear any previous NER data to ensure a clean run
        with driver.session() as session:
            logging.info("Clearing any previous NER-generated data...")
            session.run("MATCH (n:Paper) DETACH DELETE n").consume()
            # This will now remove all nodes that were created *only* by the NER script.
            session.run("MATCH (e:Entity) WHERE size(labels(e)) = 1 DETACH DELETE e").consume()


        extract_and_load_entities(driver)

        driver.close()

    except Exception as e:
        logging.error(f"An error occurred in the main pipeline: {e}")

if __name__ == "__main__":
    main()
