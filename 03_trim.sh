#!/bin/bash
# =============================================================================
# 03_trim.sh
# Trim adapters and low-quality bases using Trimmomatic (paired-end mode)
# Adapters: TruSeq3-PE-2 + NEB kit (TruSeq3-PE-2_NEB.fa)
# =============================================================================

set -euo pipefail

TRIMMOMATIC="/opt/local/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar"
ADAPTERS="TruSeq3-PE-2_NEB.fa"
RAW_DIR="SRA_all_raw_fastq"
TRIM_DIR="trimmed_reads"
THREADS=16

mkdir -p "$TRIM_DIR"

echo "=== Trimming paired-end reads with Trimmomatic ==="

for infile in "${RAW_DIR}"/*_1.fastq; do
    base=$(basename "${infile}" _1.fastq)
    echo "Processing: $base"

    java -jar "$TRIMMOMATIC" PE -threads "$THREADS" \
        "${RAW_DIR}/${base}_1.fastq" \
        "${RAW_DIR}/${base}_2.fastq" \
        "${TRIM_DIR}/${base}_R1_PE.fastq" \
        "${TRIM_DIR}/${base}_R1_SE.fastq" \
        "${TRIM_DIR}/${base}_R2_PE.fastq" \
        "${TRIM_DIR}/${base}_R2_SE.fastq" \
        ILLUMINACLIP:"${ADAPTERS}":2:30:10 \
        LEADING:3 \
        TRAILING:3 \
        MAXINFO:40:0.4 \
        MINLEN:50
done

echo "=== Trimming complete. Paired-end output in: $TRIM_DIR ==="
# SE (single-end/unpaired) files are kept but not used downstream
