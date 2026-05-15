#!/bin/bash
# =============================================================================
# 05_htseq_count.sh
# Count reads per gene using HTSeq-count on name-sorted BAM files
# Annotation: GCF_000002035.6 genomic GTF (GRCz11)
# =============================================================================

set -euo pipefail

SAMTOOLS="/usr/local/samtools/current/samtools"
GTF="ncbi_dataset/data/GCF_000002035.6/genomic.gtf"
BAM_DIR="bam_hisat_mapped"
NAME_SORTED_DIR="bam_hisat_name_sorted"
COUNTS_DIR="htseq_counts"
THREADS=8

mkdir -p "$NAME_SORTED_DIR" "$COUNTS_DIR"

# --- Step 1: Re-sort BAM files by read name (required for HTSeq) ---
echo "=== Name-sorting BAM files ==="
for infile in "${BAM_DIR}"/*_mapped.sam; do
    base=$(basename "${infile}" _mapped.sam)
    echo "Name-sorting: $base"

    "$SAMTOOLS" view -@"$THREADS" -Sb "${BAM_DIR}/${base}_mapped.sam" \
        | "$SAMTOOLS" sort -@"$THREADS" -O bam -n \
        -o "${NAME_SORTED_DIR}/${base}_name_sorted.bam"
done

# --- Step 2: Count reads with HTSeq ---
echo "=== Running HTSeq-count ==="
for infile in "${NAME_SORTED_DIR}"/*_name_sorted.bam; do
    base=$(basename "${infile}" _name_sorted.bam)
    echo "Counting: $base"

    htseq-count \
        -f bam \
        -r name \
        "$infile" \
        "$GTF" \
        > "${COUNTS_DIR}/${base}_htseq_counts.txt"
done

echo "=== HTSeq counting complete. Count files in: $COUNTS_DIR ==="
# Output: one 2-column tab-separated file per sample (gene_id | count)
