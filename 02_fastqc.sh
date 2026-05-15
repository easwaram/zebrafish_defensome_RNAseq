#!/bin/bash
# =============================================================================
# 02_fastqc.sh
# Run FastQC on raw and trimmed reads
# =============================================================================

set -euo pipefail

THREADS=8
RAW_DIR="SRA_all_raw_fastq"
TRIMMED_DIR="trimmed_reads"
RAW_QC_OUT="fastqc_all_raw_reads"
TRIMMED_QC_OUT="fastqc_trimmed_paired_reads"

mkdir -p "$RAW_QC_OUT" "$TRIMMED_QC_OUT"

echo "=== FastQC on raw reads ==="
fastqc -t "$THREADS" "${RAW_DIR}"/*.fastq -o "$RAW_QC_OUT"

echo "=== FastQC on trimmed paired-end reads ==="
# Run after trimmomatic (03_trim.sh); only QC the PE output files
fastqc -t "$THREADS" "${TRIMMED_DIR}"/*_PE.fastq -o "$TRIMMED_QC_OUT"

echo "=== FastQC complete. Review MultiQC report for summary. ==="
# Optional: run multiqc for aggregated report
# multiqc "$RAW_QC_OUT" -o "${RAW_QC_OUT}/multiqc"
# multiqc "$TRIMMED_QC_OUT" -o "${TRIMMED_QC_OUT}/multiqc"
