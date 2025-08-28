#!/bin/bash
# download_reference_grch38.sh: Download GRCh38 genome, transcriptome, and annotation files for pipeline
# Usage: bash scripts/download_reference_grch38.sh

set -euo pipefail

REF_DIR="data/reference"
mkdir -p "$REF_DIR"

# Download GRCh38 primary assembly (GENCODE)
GENOME_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/GRCh38.primary_assembly.genome.fa.gz"
TRANSCRIPTS_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.transcripts.fa.gz"
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz"

curl --output "$REF_DIR/genome.fa.gz" "$GENOME_URL"
curl --output "$REF_DIR/transcripts.fa.gz" "$TRANSCRIPTS_URL"
curl --output "$REF_DIR/reference.gtf.gz" "$GTF_URL"

# Unzip files
gunzip -f "$REF_DIR/genome.fa.gz"
gunzip -f "$REF_DIR/transcripts.fa.gz"
gunzip -f "$REF_DIR/reference.gtf.gz"

echo "Reference files downloaded to $REF_DIR:"
echo "  - genome.fa"
echo "  - transcripts.fa"
echo "  - reference.gtf"
