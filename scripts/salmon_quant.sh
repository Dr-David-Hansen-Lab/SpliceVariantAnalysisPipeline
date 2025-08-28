#!/bin/bash
# salmon_quant.sh: Run Salmon quantification for all samples listed in config/samples.csv
set -euo pipefail

CONFIG=config/config.yaml
SAMPLES=config/samples.csv
SALMON_INDEX=$(yq '.salmon_index' $CONFIG)
THREADS=$(yq '.salmon.threads' $CONFIG)
TRIMMED_DIR=$(yq '.trimmed_dir' $CONFIG)
SALMON_QUANT_DIR=$(yq '.salmon_quant_dir' $CONFIG)

mkdir -p "$SALMON_QUANT_DIR"

while IFS=, read -r sample condition fastq1 fastq2; do
    [[ "$sample" == "sample" ]] && continue  # skip header
    outdir="$SALMON_QUANT_DIR/$sample"
    mkdir -p "$outdir"
    salmon quant -i "$SALMON_INDEX" -l A -1 "$TRIMMED_DIR/${sample}_R1.trimmed.fastq.gz" -2 "$TRIMMED_DIR/${sample}_R2.trimmed.fastq.gz" -p "$THREADS" -o "$outdir"
done < "$SAMPLES"
