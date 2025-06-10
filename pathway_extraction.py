import pandas as pd
from neo4j import GraphDatabase
import logging
import os
import urllib.request
import pinecone
from sentence_transformers import SentenceTransformer
from Bio import Entrez
import time

# --- Configuration ---
# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Neo4j Configuration
NEO4J_URI = os.getenv("NEO4J_URI", "bolt://localhost:7687")
NEO4J_USER = os.getenv("NEO4J_USER", "neo4j")
NEO4J_PASSWORD = os.getenv("NEO4J_PASSWORD", "your_password") # Change this

# Pinecone Configuration
PINECONE_API_KEY = os.getenv("PINECONE_API_KEY", "YOUR_PINECONE_API_KEY")
PINECONE_ENVIRONMENT = os.getenv("PINECONE_ENVIRONMENT", "gcp-starter") # Or your environment
PINECONE_INDEX_NAME = "gwas-rag-index"

# PubMed Configuration
Entrez.email = "your.email@example.com"  # Be polite to NCBI servers

# --- Data Download Functions ---
def download_file(url, directory, filename):
    """Downloads a file from a URL if it doesn't already exist."""
    if not os.path.exists(directory):
        os.makedirs(directory)
    filepath = os.path.join(directory, filename)
    if not os.path.exists(filepath):
        logging.info(f"Downloading {filename} from {url}...")
        urllib.request.urlretrieve(url, filepath)
        logging.info(f"Downloaded {filename} successfully.")
    else:
        logging.info(f"{filename} already exists. Skipping download.")
    return filepath

# --- Neo4j Data Loading Functions ---

def import_ontology(driver, local_path, file_format='RDF/XML'):
    """Imports an ontology (OWL/OBO) into Neo4j using Neosemantics."""
    try:
        with driver.session() as session:
            # First, set up graph configurations for neosemantics
            session.run("CALL n10s.graphconfig.init({ handleVocabUris: 'IGNORE' });")
            logging.info(f"Importing {local_path} into Neo4j...")
            # Use file:/// protocol for local files
            query = f"CALL n10s.rdf.import.fetch('file:///{os.path.abspath(local_path)}', '{file_format}')"
            session.run(query)
            logging.info(f"Successfully imported {local_path}.")
    except Exception as e:
        logging.error(f"Failed to import ontology from {local_path}: {e}")

def import_reactome_pathways(driver, filepath):
    """Parses Reactome pathway file and links them to genes."""
    try:
        df = pd.read_csv(filepath, sep='\t', header=None, names=['UniProt_ID', 'Reactome_ID', 'URL', 'Pathway_Name', 'Evidence', 'Species'])
        df = df[df['Species'] == 'Homo sapiens'] # Filter for human pathways

        with driver.session() as session:
            logging.info("Importing Reactome pathways...")
            for _, row in df.iterrows():
                query = """
                MATCH (g:Gene {id: $uniprot_id})
                MERGE (p:Pathway:Reactome {id: $pathway_id})
                ON CREATE SET p.name = $pathway_name, p.url = $url
                MERGE (g)-[:PARTICIPATES_IN]->(p)
                """
                session.run(query, uniprot_id=row['UniProt_ID'], pathway_id=row['Reactome_ID'], pathway_name=row['Pathway_Name'], url=row['URL'])
            logging.info("Reactome pathways imported successfully.")
    except Exception as e:
        logging.error(f"Error importing Reactome pathways: {e}")

def fetch_and_load_pubmed_articles(driver, gene_symbols, max_articles=5):
    """Fetches PubMed articles for a list of genes and loads them into Neo4j."""
    logging.info(f"Fetching PubMed articles for {len(gene_symbols)} genes...")
    for symbol in gene_symbols:
        try:
            # Search PubMed
            handle = Entrez.esearch(db="pubmed", term=f"{symbol}[Gene Name] AND Homo sapiens[Organism]", retmax=max_articles)
            record = Entrez.read(handle)
            handle.close()
            pmids = record["IdList"]

            if not pmids:
                continue

            # Fetch article details
            handle = Entrez.efetch(db="pubmed", id=pmids, rettype="medline", retmode="xml")
            articles = Entrez.read(handle)["PubmedArticle"]
            handle.close()

            with driver.session() as session:
                for article in articles:
                    try:
                        article_info = article['MedlineCitation']['Article']
                        pmid = article['MedlineCitation']['PMID']
                        title = article_info['ArticleTitle']
                        abstract = "".join(article_info.get('Abstract', {}).get('AbstractText', []))

                        if not abstract: # Skip articles without abstracts
                            continue

                        query = """
                        MATCH (g:Gene {symbol: $gene_symbol})
                        MERGE (a:Article:PubMed {id: $pmid})
                        ON CREATE SET a.title = $title, a.abstract = $abstract
                        MERGE (g)-[:MENTIONED_IN]->(a)
                        """
                        session.run(query, gene_symbol=symbol, pmid=pmid, title=title, abstract=abstract)
                    except Exception as e:
                        logging.warning(f"Could not process article for gene {symbol}. Reason: {e}")
            time.sleep(1) # Be courteous to NCBI APIs
        except Exception as e:
            logging.error(f"Failed to fetch articles for gene {symbol}: {e}")
    logging.info("Finished fetching PubMed articles.")


# --- Pinecone Integration Functions ---

def initialize_pinecone():
    """Initializes and returns a Pinecone index."""
    pinecone.init(api_key=PINECONE_API_KEY, environment=PINECONE_ENVIRONMENT)
    if PINECONE_INDEX_NAME not in pinecone.list_indexes():
        logging.info(f"Creating Pinecone index: {PINECONE_INDEX_NAME}")
        # The embedding model dimension for 'all-MiniLM-L6-v2' is 384
        pinecone.create_index(PINECONE_INDEX_NAME, dimension=384, metric="cosine")
    return pinecone.Index(PINECONE_INDEX_NAME)

def vectorize_and_upsert(index, model, driver):
    """Fetches data from Neo4j, creates embeddings, and upserts to Pinecone."""
    with driver.session() as session:
        # Get data from Neo4j
        result = session.run("""
        MATCH (n) WHERE (n:GO_Term OR n:Pathway OR n:Article) AND n.name IS NOT NULL
        RETURN elementId(n) as id, labels(n)[0] as type, COALESCE(n.name, n.title) as text, n.abstract as abstract
        """)
        records = list(result) # Consume the result fully

    logging.info(f"Vectorizing {len(records)} nodes from Neo4j...")
    vectors_to_upsert = []
    for record in records:
        # Create a descriptive text for embedding
        content = record['text']
        if record['abstract']:
            content += ". " + record['abstract']

        embedding = model.encode(content).tolist()
        metadata = {"type": record['type'], "text": record['text']}
        vectors_to_upsert.append((record['id'], embedding, metadata))

        # Upsert in batches to Pinecone
        if len(vectors_to_upsert) >= 100:
            index.upsert(vectors=vectors_to_upsert)
            vectors_to_upsert = []

    # Upsert any remaining vectors
    if vectors_to_upsert:
        index.upsert(vectors=vectors_to_upsert)

    logging.info("Finished upserting vectors to Pinecone.")


# --- Main Execution ---
def main():
    # --- 1. Download External Data ---
    logging.info("--- Step 1: Downloading data ---")
    download_dir = "downloads"
    # Gene Ontology
    go_owl_url = "http://purl.obolibrary.org/obo/go/go-basic.owl"
    go_owl_path = download_file(go_owl_url, download_dir, "go-basic.owl")
    # Human Phenotype Ontology
    hpo_owl_url = "http://purl.obolibrary.org/obo/hp.owl"
    hpo_owl_path = download_file(hpo_owl_url, download_dir, "hp.owl")
    # Reactome Pathways
    reactome_url = "https://reactome.org/download/current/UniProt2Reactome.txt"
    reactome_path = download_file(reactome_url, download_dir, "UniProt2Reactome.txt")

    # --- 2. Initialize Connections ---
    logging.info("--- Step 2: Initializing connections ---")
    neo4j_driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASSWORD))
    pinecone_index = initialize_pinecone()
    embedding_model = SentenceTransformer('all-MiniLM-L6-v2')


    # --- 3. Load Data into Neo4j ---
    logging.info("--- Step 3: Loading data into Neo4j ---")
    # Import ontologies
    import_ontology(neo4j_driver, go_owl_path)
    import_ontology(neo4j_driver, hpo_owl_path)
    # Import pathways (your existing GAF/VCF import functions would go here too)
    # ... (call your existing functions for GAF, VCF)
    import_reactome_pathways(neo4j_driver, reactome_path)

    # Fetch PubMed data for genes in the graph
    with neo4j_driver.session() as session:
        gene_result = session.run("MATCH (g:Gene) WHERE g.symbol IS NOT NULL RETURN collect(g.symbol) as symbols")
        gene_symbols = gene_result.single()['symbols']
    fetch_and_load_pubmed_articles(neo4j_driver, gene_symbols)

    # --- 4. Vectorize and Store in Pinecone ---
    logging.info("--- Step 4: Vectorizing data for Pinecone ---")
    vectorize_and_upsert(pinecone_index, embedding_model, neo4j_driver)


    # --- 5. Clean up ---
    logging.info("--- Step 5: Cleaning up ---")
    neo4j_driver.close()
    logging.info("Pipeline finished successfully!")

if __name__ == "__main__":
    main()
