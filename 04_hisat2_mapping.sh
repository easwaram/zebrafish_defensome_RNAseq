#!/bin/bash
# =============================================================================
# 04_hisat2_mapping.sh
# Splice-aware alignment to Danio rerio genome (GRCz11) using HISAT2
# Reference genome: GCF_000002035.6 (downloaded from NCBI)
# =============================================================================

set -euo pipefail

HISAT2="/usr/local/hisat/hisat2"
HISAT2_BUILD="/usr/local/hisat/hisat2-build"
SAMTOOLS="/usr/local/samtools/current/samtools"

GENOME_FA="ncbi_dataset/data/GCF_000002035.6/GCF_000002035.6_GRCz11_genomic.fna"
INDEX_DIR="hisat2_index"
INDEX_PREFIX="${INDEX_DIR}/drerio_index_GRCz11"
TRIMMED_DIR="trimmed_reads"
SAM_DIR="hisat_mapped_reads"
BAM_DIR="bam_hisat_mapped"
THREADS=8

mkdir -p "$INDEX_DIR" "$SAM_DIR" "$BAM_DIR"

# --- Step 1: Download reference genome ---
echo "=== Downloading GRCz11 reference genome ==="
datasets download genome accession GCF_000002035.6
unzip ncbi_dataset.zip

# --- Step 2: Build HISAT2 index ---
echo "=== Building HISAT2 index ==="
"$HISAT2_BUILD" -p 6 -f "$GENOME_FA" "$INDEX_PREFIX"

# --- Step 3: Align trimmed paired-end reads ---
echo "=== Aligning reads with HISAT2 ==="
for infile in "${TRIMMED_DIR}"/*_R1_PE.fastq; do
    base=$(basename "${infile}" _R1_PE.fastq)
    echo "Mapping: $base"

    "$HISAT2" -x "$INDEX_PREFIX" \
        -1 "${TRIMMED_DIR}/${base}_R1_PE.fastq" \
        -2 "${TRIMMED_DIR}/${base}_R2_PE.fastq" \
        -p "$THREADS" \
        -S "${SAM_DIR}/${base}_mapped.sam"
done

# --- Step 4: Convert SAM to sorted BAM using SAMtools ---
echo "=== Converting SAM to sorted BAM ==="
for infile in "${SAM_DIR}"/*_mapped.sam; do
    base=$(basename "${infile}" _mapped.sam)
    echo "Converting: $base"

    "$SAMTOOLS" view -@"$THREADS" -Sb "${SAM_DIR}/${base}_mapped.sam" \
        | "$SAMTOOLS" sort -@"$THREADS" -O bam \
        -o "${BAM_DIR}/${base}_sorted.bam"
done

echo "=== Alignment and BAM sorting complete. Summary logs in: hisat_summary_alignment/ ==="
