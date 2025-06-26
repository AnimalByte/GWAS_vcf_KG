import os
import logging
import re
from neo4j import GraphDatabase
import chromadb
from sentence_transformers import SentenceTransformer
from llama_cpp import Llama
import spacy

# --- Configuration ---
# Sets up logging to monitor the script's progress.
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# --- Paths and Model Configuration ---
# Neo4j connection details
NEO4J_URI = os.getenv("NEO4J_URI", "bolt://localhost:7687")
NEO4J_USER = os.getenv("NEO4J_USER", "neo4j")
NEO4J_PASSWORD = os.getenv("NEO4J_PASSWORD", "password")

# ChromaDB connection details
DB_PATH = "chroma_db"
COLLECTION_NAME = "pubmed_abstracts"

# Embedding Model (must be the same one used to create the database)
EMBEDDING_MODEL_NAME = 'pritamdeka/S-PubMedBert-MS-MARCO'

# NER Model for entity extraction
SPACY_MODEL_NAME = "en_core_sci_lg"

# Local LLM Configuration
MODEL_PATH = "./models/Phi-3-medium-128k-instruct-Q6_K_L.gguf" # Using the highest quality model

# --- Hardware Configuration ---
N_GPU_LAYERS = -1
N_CTX = 4096

def flatten_and_unique(items_list):
    """
    Helper function to flatten a list that may contain sublists and return
    a unique, case-insensitively filtered list of strings.
    """
    if not items_list:
        return []
    
    flat_list = []
    for item in items_list:
        if isinstance(item, list):
            flat_list.extend(item)
        elif item is not None:
            flat_list.append(item)
    
    seen = set()
    unique_list = []
    for item in flat_list:
        if str(item).lower() not in seen:
            seen.add(str(item).lower())
            unique_list.append(str(item))
    return unique_list


class GraphRAGQueryEngine:
    def __init__(self):
        logging.info("Initializing GraphRAG Query Engine...")
        self.neo4j_driver = GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASSWORD))
        
        self.embedding_model = SentenceTransformer(EMBEDDING_MODEL_NAME)
        self.chroma_client = chromadb.PersistentClient(path=DB_PATH)
        self.collection = self.chroma_client.get_collection(name=COLLECTION_NAME)
        
        # Load the scispaCy NER model
        logging.info(f"Loading NER model: {SPACY_MODEL_NAME}...")
        try:
            spacy.prefer_gpu()
            self.nlp = spacy.load(SPACY_MODEL_NAME)
            logging.info("NER model loaded successfully on GPU (if available).")
        except OSError:
            logging.error(f"Model '{SPACY_MODEL_NAME}' not found. Please run the correct installation commands.")
            exit()
            
        logging.info(f"Loading LLM from path: {MODEL_PATH}")
        if not os.path.exists(MODEL_PATH):
            logging.error(f"LLM model not found at {MODEL_PATH}.")
            exit()
            
        self.llm = Llama(model_path=MODEL_PATH, n_gpu_layers=N_GPU_LAYERS, n_ctx=N_CTX, verbose=False)
        logging.info("LLM loaded successfully.")

        # Pre-fetch graph entity names for robust matching
        self._build_entity_maps()

    def _build_entity_maps(self):
        """Fetches all gene, phenotype, and disease names from the graph to create an in-memory map."""
        logging.info("Building in-memory maps for entity resolution...")
        self.phenotype_map = {}
        self.gene_map = {}
        self.pathway_map = {}
        try:
            with self.neo4j_driver.session() as session:
                # Fetch Phenotypes and Diseases
                pheno_result = session.run("MATCH (p) WHERE (p:Phenotype OR p:Disease) AND p.name IS NOT NULL RETURN p.name, labels(p)[0] AS label")
                for record in pheno_result:
                    name_field = record["p.name"]
                    names = name_field if isinstance(name_field, list) else [name_field]
                    for name in names:
                        self.phenotype_map[name.lower()] = {"name": name, "type": record["label"]}
                
                # Fetch Genes
                gene_result = session.run("MATCH (g:Gene) WHERE g.symbol IS NOT NULL RETURN g.symbol")
                for record in gene_result:
                    symbol = record["g.symbol"]
                    self.gene_map[symbol.lower()] = symbol

                # Fetch Pathways
                pathway_result = session.run("MATCH (p:Pathway) WHERE p.name IS NOT NULL RETURN p.name")
                for record in pathway_result:
                    pathway_name = record["p.name"]
                    self.pathway_map[pathway_name.lower()] = pathway_name

            logging.info(f"Built maps with {len(self.phenotype_map)} phenotypes/diseases, {len(self.gene_map)} genes, and {len(self.pathway_map)} pathways.")
        except Exception as e:
            logging.error(f"Failed to build entity maps: {e}")


    def close(self):
        self.neo4j_driver.close()

    def retrieve_context(self, question):
        """
        The main retrieval function. It uses NER to identify entities, resolves them against
        the graph, queries for context, and then performs an expanded vector search.
        """
        doc = self.nlp(question)
        entities = [{"text": ent.text.lower(), "label": ent.label_} for ent in doc.ents]
        logging.info(f"NER found potential entities: {entities}")

        graph_context = "No relevant data found in the graph."
        full_context = ""
        found_specific_entity = False

        if entities:
            with self.neo4j_driver.session() as session:
                # **UPGRADED**: Loop through all found entities and build a context for each resolved one.
                for entity in entities:
                    entity_text = entity['text']
                    
                    # --- Resolve against our in-memory maps ---
                    if entity_text in self.gene_map:
                        found_specific_entity = True
                        gene_symbol = self.gene_map[entity_text]
                        logging.info(f"Resolved '{entity_text}' as Gene: {gene_symbol}")
                        context_query = """
                        MATCH (g:Gene {symbol: $symbol})
                        OPTIONAL MATCH (g)-[:PARTICIPATES_IN]->(p:Pathway)
                        OPTIONAL MATCH (g)-[:ASSOCIATED_WITH]->(h:Phenotype)
                        OPTIONAL MATCH (g)<-[:AFFECTS]-(m:Mutation)
                        RETURN g.symbol AS name, 'Gene' AS type, 
                               collect(DISTINCT p.name) as pathways,
                               collect(DISTINCT h.name) as phenotypes,
                               collect(DISTINCT {id: m.id, significance: m.clinical_significance})[..5] as mutations
                        """
                        result = session.run(context_query, symbol=gene_symbol).single()
                        if result:
                            full_context += f"\nEntity: {result['name']} (Type: {result['type']})\n"
                            if result.get('pathways'): full_context += f"  - Associated Pathways: {', '.join(flatten_and_unique(result['pathways'])[:3])}\n"
                            if result.get('phenotypes'): full_context += f"  - Associated Phenotypes: {', '.join(flatten_and_unique(result['phenotypes'])[:3])}\n"

                    elif entity_text in self.phenotype_map:
                        found_specific_entity = True
                        phenotype_data = self.phenotype_map[entity_text]
                        phenotype_name = phenotype_data['name']
                        phenotype_type = phenotype_data['type']
                        logging.info(f"Resolved '{entity_text}' as {phenotype_type}: {phenotype_name}")
                        context_query = "MATCH (p {name: $name})<-[*]-(g:Gene) RETURN collect(DISTINCT g.symbol)[..10] as genes"
                        result = session.run(context_query, name=phenotype_name).single()
                        if result:
                            full_context += f"Entity: {phenotype_name} (Type: {phenotype_type})\n"
                            if result.get('genes'): full_context += f"  - Associated Genes: {', '.join(result['genes'])}\n"
                    
                    elif entity_text in self.pathway_map:
                        found_specific_entity = True
                        pathway_name = self.pathway_map[entity_text]
                        logging.info(f"Resolved '{entity_text}' as Pathway: {pathway_name}")
                        context_query = "MATCH (p:Pathway {name: $name})<-[:PARTICIPATES_IN]-(g:Gene) RETURN collect(DISTINCT g.symbol)[..10] as genes"
                        result = session.run(context_query, name=pathway_name).single()
                        if result:
                            full_context += f"Entity: {pathway_name} (Type: Pathway)\n"
                            if result.get('genes'): full_context += f"  - Associated Genes: {', '.join(result['genes'])}\n"


        # If NER finds no specific, resolvable entities, fall back to the discovery query
        if not found_specific_entity:
            logging.info("No specific entities resolved. Running discovery query...")
            with self.neo4j_driver.session() as session:
                discovery_query = """
                MATCH (m:Mutation)-[:AFFECTS]->(g:Gene)
                WHERE m.impact = 'HIGH' OR m.impact = 'MODERATE'
                WITH g, count(m) AS mutation_count
                ORDER BY mutation_count DESC LIMIT 5
                OPTIONAL MATCH (g)-[:PARTICIPATES_IN]->(p:Pathway)
                OPTIONAL MATCH (g)-[:ASSOCIATED_WITH]->(h:Phenotype)
                RETURN g.symbol AS gene, 
                       collect(DISTINCT p.name) AS pathways, 
                       collect(DISTINCT h.name) AS phenotypes
                """
                result = session.run(discovery_query)
                records = list(result)
                if records:
                    full_context = "Top impactful genes from the dataset and their key associations:\n"
                    for record in records:
                        full_context += f"- Gene: {record['gene']}\n"
                        if record['pathways']: full_context += f"  - Pathways: {', '.join(flatten_and_unique(record['pathways'])[:3])}\n"
                        if record['phenotypes']: full_context += f"  - Phenotypes: {', '.join(flatten_and_unique(record['phenotypes'])[:3])}\n"

        if full_context:
            graph_context = full_context.strip()
        
        expanded_query = question + " " + graph_context
        logging.info(f"Expanded vector search query: '{expanded_query}'")
        vector_context = self.retrieve_from_vector(expanded_query)
        
        return graph_context, vector_context

    def retrieve_from_vector(self, query_text):
        """Retrieves relevant abstracts from ChromaDB based on semantic similarity."""
        logging.info("Querying vector database for relevant abstracts...")
        query_embedding = self.embedding_model.encode(query_text)
        results = self.collection.query(query_embeddings=[query_embedding.tolist()], n_results=3)
        context = ""
        for i, doc in enumerate(results['documents'][0]):
            context += f"--- Abstract {i+1} (PMID: {results['ids'][0][i]}) ---\n{doc}\n\n"
        return context

    def answer_question(self, question):
        """Main RAG pipeline to answer a user's question."""
        logging.info(f"Received question: {question}")
        graph_context, vector_context = self.retrieve_context(question)
        prompt = f"""
        [INST] You are a helpful biomedical research assistant.
        Your task is to answer the following question based *only* on the context provided.
        Synthesize the information from the structured 'Graph Context' and the unstructured 'Vector Context' from PubMed abstracts.
        Be concise and cite the PubMed IDs (PMID) of any abstracts you use in your answer.

        **Graph Context (Structured Data):**
        {graph_context}

        **Vector Context (Unstructured Abstracts):**
        {vector_context}
        
        **Question:**
        {question} [/INST]
        """
        logging.info("Generating answer with LLM...")
        response = self.llm(prompt, max_tokens=512, stream=True)
        
        print("--- Answer ---")
        for chunk in response:
            print(chunk['choices'][0]['text'], end="", flush=True)
        print("\n")

if __name__ == "__main__":
    engine = GraphRAGQueryEngine()
    
    print("\n--- Bio-KG RAG System (Final Engine) ---")
    print("Ask a question about a gene, phenotype, or pathway (e.g., 'What genes are associated with Alopecia?') or type 'exit' to quit.")
    while True:
        user_question = input("> ")
        if user_question.lower() == 'exit':
            break
        engine.answer_question(user_question)

    engine.close()
    print("Session ended.")
