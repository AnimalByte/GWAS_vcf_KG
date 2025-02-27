import pandas as pd
from neo4j import GraphDatabase
import logging
import os

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# -------------------- 🧠 Import OWL Ontology into Neo4j --------------------
def import_owl_to_neo4j(driver, owl_file_path):
    """
    Import an OWL ontology into Neo4j using Neosemantics (n10s).
    """
    if not os.path.exists(owl_file_path):
        logging.error(f"OWL file not found at path: {owl_file_path}")
        return

    try:
        with driver.session() as session:
            logging.info("Importing OWL file into Neo4j using n10s...")
            query = """
            CALL n10s.rdf.import.fetch($owl_file, 'RDF/XML', {})
            """
            session.run(query, owl_file=owl_file_path)
            logging.info("OWL file imported successfully.")
    except Exception as e:
        logging.error(f"Failed to import OWL file: {e}")

# -------------------- 🔬 Parse GAF File --------------------
def parse_gaf_file(gaf_file_path):
    """
    Parse a GAF file with 17 columns (GAF 2.2 specification).
    Keeps all columns for downstream tasks.
    """
    if not os.path.exists(gaf_file_path):
        logging.error(f"GAF file not found at path: {gaf_file_path}")
        return None

    try:
        logging.info(f"Parsing GAF file: {gaf_file_path}")

        # Read the GAF file, ignoring comment lines (lines starting with "!")
        df = pd.read_csv(
            gaf_file_path,
            sep="\t",
            comment="!",
            header=None,
            compression="gzip",
            low_memory=False
        )

        # Define column names for GAF 2.2 format (17 columns)
        gaf_columns = [
            "DB", "DB_Object_ID", "DB_Object_Symbol", "Relation", "GO_ID",
            "DB_Reference", "Evidence_Code", "With_From", "Aspect",
            "DB_Object_Name", "DB_Object_Synonym", "DB_Object_Type",
            "Taxon", "Date", "Assigned_By", "Annotation_Extension", "Gene_Product_Form_ID"
        ]

        # Check if column count matches GAF 2.2 format
        if df.shape[1] != 17:
            logging.error(f"Unexpected column count ({df.shape[1]} instead of 17). Check file structure!")
            return None

        # Assign column names
        df.columns = gaf_columns

        logging.info(f"Successfully parsed {len(df)} rows from GAF file.")
        return df

    except Exception as e:
        logging.error(f"Error reading GAF file: {e}")
        return None

# -------------------- 🏷️ Import Gene-GO Annotations into Neo4j --------------------
def import_annotations_to_neo4j(driver, annotations):
    """
    Import gene-GO term annotations into Neo4j, ensuring Gene symbols are assigned.
    """
    try:
        def import_annotation(tx, annotation):
            query = """
            MERGE (g:Gene {id: $DB_Object_ID})
            SET g.symbol = $DB_Object_Symbol
            MERGE (go:GO_Term {id: $GO_ID})
            MERGE (g)-[:ANNOTATED_WITH {evidence: $Evidence_Code}]->(go)
            """
            tx.run(query, **annotation)

        with driver.session() as session:
            logging.info("Importing GAF annotations into Neo4j...")
            for _, annotation in annotations.iterrows():
                session.execute_write(import_annotation, annotation.to_dict())
            logging.info("Annotations imported successfully.")
    except Exception as e:
        logging.error(f"An error occurred while importing annotations: {e}")

# -------------------- 🔬 Import VCF Mutations into Neo4j --------------------
def import_vcf_mutations(driver, vcf_annotation_path):
    """
    Parse the VEP annotated VCF file and import mutation data into Neo4j.
    Each mutation is connected to its gene using the 'Gene' field.
    """
    try:
        # Define the expected column names (since the header line is commented out)
        vcf_columns = [
            "Uploaded_variation", "Location", "Allele", "Gene",
            "Feature", "Feature_type", "Consequence", "SIFT", "PolyPhen", "IMPACT"
        ]

        # Read the VEP annotated file with the defined columns.
        df = pd.read_csv(vcf_annotation_path, sep="\t", comment="#", names=vcf_columns)
        if df.empty:
            logging.error("The VEP annotated VCF file is empty!")
            return

        with driver.session() as session:
            logging.info("Importing VCF mutations into Neo4j...")
            for _, row in df.iterrows():
                mutation = row.to_dict()
                # Optionally, strip whitespace from the keys (if necessary)
                mutation = {k.strip(): v for k, v in mutation.items()}
                query = """
                MERGE (m:Mutation {id: $Uploaded_variation})
                SET m.location = $Location,
                    m.allele = $Allele,
                    m.feature = $Feature,
                    m.feature_type = $Feature_type,
                    m.consequence = $Consequence,
                    m.impact = $IMPACT,
                    m.sift = $SIFT,
                    m.polyphen = $PolyPhen
                WITH m
                MATCH (g:Gene)
                WHERE g.symbol = $Gene
                MERGE (m)-[:AFFECTS]->(g)
                """
                session.run(query, **mutation)
            logging.info("VCF mutations imported successfully.")
    except Exception as e:
        logging.error(f"Error importing VCF mutations: {e}")

# -------------------- 🚀 Main Execution --------------------
def main():
    try:
        # File paths
        owl_file_path = "downloads/go-plus.owl"
        gaf_file_path = "downloads/goa_human.gaf.gz"
        vcf_annotation_path = "results/vep_annotated.txt"

        # Establish Neo4j connection (update authentication if needed)
        uri = "bolt://localhost:7687"
        driver = GraphDatabase.driver(uri, auth=None)

        # 1. Import GO Ontology from the OWL file
        import_owl_to_neo4j(driver, owl_file_path)

        # 2. Import Gene-GO Annotations from the GAF file
        annotations = parse_gaf_file(gaf_file_path)
        if annotations is not None and not annotations.empty:
            import_annotations_to_neo4j(driver, annotations)
        else:
            logging.error("No annotations found in the GAF file or the file is empty.")

        # 3. Import VCF Mutations from the annotated VEP file
        import_vcf_mutations(driver, vcf_annotation_path)

        # Close Neo4j connection
        driver.close()

    except Exception as e:
        logging.error(f"An error occurred: {e}")

if __name__ == "__main__":
    main()
