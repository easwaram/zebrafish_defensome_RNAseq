#!/bin/bash
# =============================================================================
# 06_salmon_quant.sh
# Pseudo-alignment and quantification using Salmon (selective alignment mode)
# Reference: GRCz11 transcriptome + genome decoy
# =============================================================================

set -euo pipefail

GENOME_FA="ncbi_dataset/data/GCF_000002035.6/GCF_000002035.6_GRCz11_genomic.fna.gz"
TRANSCRIPTOME_FA="ncbi_dataset/data/GCF_000002035.6/rna.fna"
SALMON_INDEX="salmon_index_GRCz11"
TRIMMED_DIR="trimmed_reads"
OUT_DIR="salmon_counts"
THREADS=16

mkdir -p "$SALMON_INDEX" "$OUT_DIR"

# --- Step 1: Build decoy-aware Salmon index ---
echo "=== Building Salmon index with genome decoy ==="

# Extract genome sequence names for decoy list
grep "^>" <(gunzip -c "$GENOME_FA") | cut -d " " -f 1 > decoys.txt
sed -i.bak -e 's/>//g' decoys.txt

# Concatenate transcriptome + genome for selective alignment
cat "$TRANSCRIPTOME_FA" <(gunzip -c "$GENOME_FA") > GRCz11_salmon_reference.fa

salmon index \
    -t GRCz11_salmon_reference.fa \
    -d decoys.txt \
    -p "$THREADS" \
    -i "$SALMON_INDEX"

# --- Step 2: Quantify each sample ---
echo "=== Quantifying with Salmon ==="
for infile in "${TRIMMED_DIR}"/*_R1_PE.fastq; do
    base=$(basename "${infile}" _R1_PE.fastq)
    echo "Quantifying: $base"

    salmon quant \
        -i "$SALMON_INDEX" \
        -l A \
        -1 "${TRIMMED_DIR}/${base}_R1_PE.fastq" \
        -2 "${TRIMMED_DIR}/${base}_R2_PE.fastq" \
        -p "$THREADS" \
        --validateMappings \
        --rangeFactorizationBins 4 \
        --seqBias \
        --gcBias \
        -o "${OUT_DIR}/${base}_quant"
done

echo "=== Salmon quantification complete. Results in: $OUT_DIR ==="
# Library type (-l A) was auto-detected as IU (inward, unstranded)
# Mapping rates are recorded in each sample's logs/salmon_quant.log
