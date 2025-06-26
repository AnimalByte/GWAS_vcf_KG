import os
import logging
import chromadb
from sentence_transformers import SentenceTransformer
import torch

# --- Configuration ---
# Sets up logging to monitor the script's progress.
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# --- Paths and Model Configuration ---
# Directory where the downloaded abstracts are stored.
ABSTRACTS_DIR = "pubmed_abstracts"
# Path to store the local vector database.
DB_PATH = "chroma_db"
# Name of the collection within the database.
COLLECTION_NAME = "pubmed_abstracts"
# This is a version of PubMedBERT adapted for creating high-quality sentence/document embeddings.
MODEL_NAME = 'pritamdeka/S-PubMedBert-MS-MARCO'

# --- Main Script ---
def create_and_store_embeddings():
    """
    Scans the abstracts, generates embeddings using a biomedical-specific
    transformer model, and stores them in a local ChromaDB vector database.
    This version includes a check to prevent duplicate IDs from being added.
    """
    logging.info("--- Starting Embedding Creation and Storage ---")

    # --- Step 1: Initialize Model ---
    # Check if a GPU is available and set the device accordingly.
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    logging.info(f"Using device: {device}")
    
    # Load the SentenceTransformer model. This will download the model weights if not cached.
    logging.info(f"Loading sentence transformer model: {MODEL_NAME}")
    model = SentenceTransformer(MODEL_NAME, device=device)
    logging.info("Model loaded successfully.")

    # --- Step 2: Initialize Vector Database and Check for Existing Entries ---
    # Set up the ChromaDB client. 'PersistentClient' saves the DB to disk.
    client = chromadb.PersistentClient(path=DB_PATH)
    
    # Get or create the collection. A collection is like a table in a traditional DB.
    logging.info(f"Initializing ChromaDB collection: '{COLLECTION_NAME}'")
    collection = client.get_or_create_collection(name=COLLECTION_NAME)

    # **FIX**: Get all IDs that are already in the database to avoid reprocessing.
    # This makes the script idempotent (safe to rerun).
    existing_ids_results = collection.get(include=[]) # Fetches minimal data, just IDs
    processed_ids = set(existing_ids_results['ids'])
    logging.info(f"Found {len(processed_ids)} existing documents in the database. They will be skipped.")

    # --- Step 3: Process Abstracts and Create Embeddings ---
    logging.info(f"Scanning for abstracts in: '{ABSTRACTS_DIR}'")
    
    # Prepare lists to hold the data for batch insertion into ChromaDB.
    documents_batch = []
    metadatas_batch = []
    ids_batch = []
    batch_size = 50  # Process files in batches for efficiency.
    
    total_files_scanned = 0
    # Walk through the directory structure.
    for gene_symbol in os.listdir(ABSTRACTS_DIR):
        gene_dir = os.path.join(ABSTRACTS_DIR, gene_symbol)
        if os.path.isdir(gene_dir):
            for filename in os.listdir(gene_dir):
                if filename.endswith(".txt"):
                    total_files_scanned += 1
                    pmid = os.path.splitext(filename)[0] # The filename is the PubMed ID.

                    # **FIX**: Check if the PMID has already been processed in this or a previous run.
                    if pmid in processed_ids:
                        continue # Skip this file as it's a duplicate.

                    filepath = os.path.join(gene_dir, filename)
                    with open(filepath, 'r', encoding='utf-8') as f:
                        abstract_text = f.read()

                    # Add the data to our lists for the current batch.
                    documents_batch.append(abstract_text)
                    metadatas_batch.append({"gene": gene_symbol, "pmid": pmid})
                    ids_batch.append(pmid)
                    processed_ids.add(pmid) # Add to our set to prevent duplicates within the same run.
                    
                    # If the batch is full, process and store it.
                    if len(documents_batch) >= batch_size:
                        logging.info(f"Processing batch of {len(documents_batch)} new abstracts...")
                        embeddings = model.encode(documents_batch, show_progress_bar=True)
                        collection.add(
                            embeddings=embeddings,
                            documents=documents_batch,
                            metadatas=metadatas_batch,
                            ids=ids_batch
                        )
                        logging.info(f"Successfully stored batch. Total unique documents in DB: {len(processed_ids)}")
                        documents_batch, metadatas_batch, ids_batch = [], [], []

    # Process any remaining files in the last, smaller batch.
    if documents_batch:
        logging.info(f"Processing final batch of {len(documents_batch)} new abstracts...")
        embeddings = model.encode(documents_batch, show_progress_bar=True)
        collection.add(
            embeddings=embeddings,
            documents=documents_batch,
            metadatas=metadatas_batch,
            ids=ids_batch
        )

    logging.info("--- Embedding Creation and Storage Complete! ---")
    logging.info(f"Total files scanned: {total_files_scanned}")
    logging.info(f"Total unique documents in database: {collection.count()}")
    logging.info(f"Vector database is stored at: '{DB_PATH}'")


if __name__ == "__main__":
    create_and_store_embeddings()

