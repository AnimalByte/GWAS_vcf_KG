#!/bin/bash

# This script cleans the project directory, removing all generated files,
# databases, and logs, returning it to a clean state.

echo "--- Starting Project Cleanup ---"

# Stop and remove docker containers to release file locks
echo "Stopping and removing Docker containers..."
docker-compose down

# Remove generated directories
echo "Removing generated directories..."
rm -rf __pycache__/
rm -rf pubmed_abstracts/
rm -rf chroma_db/
rm -rf backups/
sudo rm -rf neo4j/data

# Remove generated files in the root directory
echo "Removing generated files..."
rm -f *.log
rm -f *.deb
rm -f *.deb.*
rm -f *.zip
rm -f *.tar.gz
rm -f *Zone.Identifier

# Remove old/temporary scripts
echo "Removing old and temporary scripts..."
rm -f parser.py
rm -f pathway_extraction.py
rm -f pathway_extraction1.py
rm -f retrival.py
rm -f enrich_genes.py
rm -f debug_uniprot_match.py
rm -f query_engine1.py
rm -f script.sh

# Re-create the essential directories
echo "Re-creating essential directories..."
mkdir -p results
mkdir -p downloads
mkdir -p neo4j/data

# Set correct permissions for the new neo4j data directory
echo "Setting permissions for Neo4j data directory..."
sudo chown -R 7474:7474 neo4j/data

echo "--- Cleanup Complete! ---"
echo "Project is now in a clean state. You can now run the full pipeline."
