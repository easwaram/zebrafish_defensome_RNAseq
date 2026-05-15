#!/bin/bash
# =============================================================================
# 01_download_sra.sh
# Download raw FASTQ files from NCBI SRA using accession lists
# Organism: Danio rerio (zebrafish) | Tissues: gut, gill
# =============================================================================

set -euo pipefail

SRATOOLKIT="/usr/local/sratoolkit/bin"
GUT_ACC="SRR_Acc_List_gut.txt"
GILL_ACC="SRR_Acc_List_gill.txt"
OUT_GUT="SRA_raw_fastq_gut"
OUT_GILL="SRA_raw_fastq_gill"

mkdir -p "$OUT_GUT" "$OUT_GILL"

echo "=== Downloading gut samples ==="
for sra_id in $(cat "$GUT_ACC"); do
    echo "Downloading $sra_id ..."
    "$SRATOOLKIT/prefetch.3.0.10" "$sra_id"
    "$SRATOOLKIT/fasterq-dump.3.0.10" --split-files "$sra_id" --outdir "$OUT_GUT"
done

echo "=== Downloading gill samples ==="
for sra_id in $(cat "$GILL_ACC"); do
    echo "Downloading $sra_id ..."
    "$SRATOOLKIT/prefetch.3.0.10" "$sra_id"
    "$SRATOOLKIT/fasterq-dump.3.0.10" --split-files "$sra_id" --outdir "$OUT_GILL"
done

echo "=== Download complete ==="
# NOTE: Files are then renamed to the convention:
# {sra_id}_{tissue}_{control#}_{read#}.fastq
# and moved to combined directory: SRA_all_raw_fastq/
