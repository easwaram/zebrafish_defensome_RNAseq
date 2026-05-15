# Zebrafish Chemical Defensome: Gut vs Gill RNA-seq Analysis

Bulk RNA-seq analysis comparing basal chemical defensome gene expression between gut and gill tissue in adult zebrafish (*Danio rerio*). This project was designed to inform the experimental design of future baits target enrichment sequencing experiments by determining whether tissue-specific differences in defensome expression exist prior to any chemical exposure.

---

## Background

The chemical defensome consists of genes encoding transcription factors, biotransformation enzymes, transporters, and antioxidant proteins that collectively mediate xenobiotic detoxification and elimination. While the liver defensome has been relatively well characterized in zebrafish, basal expression in tissues directly involved in xenobiotic uptake — such as the gut and gill — is largely unexplored.

This analysis uses publicly available control RNA-seq data from Xue et al. (2021) (NCBI BioProject: [PRJNA662254](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA662254) and [PRJNA670521](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA670521)), consisting of non-exposed adult male zebrafish sampled from gut and gill tissue. Repurposing existing public data avoided the need for a new experiment while providing suitable controls for a basal expression comparison.

---

## Biological Question

**Do gut and gill tissue show tissue-specific differences in the basal expression of chemical defensome genes in adult zebrafish?**

---

## Samples

| Tissue | Replicate | SRA Run |
|--------|-----------|---------|
| Gill | 1 | SRR12876848 |
| Gill | 2 | SRR12876847 |
| Gill | 3 | SRR12876838 |
| Gut | 1 | SRR12614074 |
| Gut | 2 | SRR12614073 |
| Gut | 3 | SRR12614064 |

All samples are pooled from 5 healthy adult male zebrafish per replicate. Library prep: NEBNext Ultra Directional RNA Library Prep Kit; sequencing: Illumina HiSeqX10.

---

## Pipeline Overview

Two parallel pipelines were implemented and compared:

**Pipeline 1 — Pseudo-mapping (Salmon + tximport + DESeq2)**
```
Raw FASTQ → FastQC → Trimmomatic → Salmon → tximport → DESeq2 → GO analysis
```

**Pipeline 2 — Splice-aware mapping (HISAT2 + HTSeq + DESeq2)**
```
Raw FASTQ → FastQC → Trimmomatic → HISAT2 → SAMtools → HTSeq-count → DESeq2 → GO analysis
```

| Script | Description |
|--------|-------------|
| `01_download_sra.sh` | Download raw FASTQ files from NCBI SRA |
| `02_fastqc.sh` | Quality control on raw and trimmed reads |
| `03_trim.sh` | Adapter and quality trimming with Trimmomatic |
| `04_hisat2_mapping.sh` | Splice-aware alignment to GRCz11 + SAM→BAM conversion |
| `05_htseq_count.sh` | Read counting per gene (Pipeline 2) |
| `06_salmon_quant.sh` | Pseudo-alignment and quantification (Pipeline 1) |
| `07_DESeq2_analysis.R` | Differential expression, LFC shrinkage, GO enrichment |

---

## Key Results

**PCA:** Clear separation of gut and gill samples, with PC1 and PC2 explaining 93% (Pipeline 1) and 97% (Pipeline 2) of variance, confirming strong tissue-specific transcriptomic differences.

**Defensome DEGs:**
- Pipeline 1 (Salmon): 314 differentially expressed defensome genes — 230 upregulated, 84 downregulated in gut vs gill
- Pipeline 2 (HISAT2): 374 differentially expressed defensome genes — 283 upregulated, 91 downregulated in gut vs gill

**GO Enrichment:**
- *Gut* upregulation: xenobiotic metabolic processes, sulfation, glutathione metabolism, cytochrome P450 activity — suggesting greater detoxification capacity in the intestine
- *Gill* upregulation: metal ion transport, protein folding, prostaglandin biosynthesis — consistent with the gill's role in ion homeostasis and environmental stress response

**Pipeline comparison:** Pipelines 1 and 2 showed strong overlap in top DEGs with highly similar log2 fold changes (e.g. *abcb4*: 8.64 vs 8.74). Pipeline 2 detected more total DEGs (12,550 vs 7,230), likely due to HISAT2's improved sensitivity for lowly abundant transcripts. Only 4 of the top 10 defensome DEGs were shared between pipelines, highlighting the importance of pipeline selection for downstream interpretation.

---

## Tools & Versions

| Tool | Version |
|------|---------|
| FastQC | 0.11.5 |
| Trimmomatic | 0.39 |
| Salmon | 1.4.0 |
| HISAT2 | 2.1.0 |
| SAMtools | 1.11 |
| HTSeq-count | 0.11.2 |
| DESeq2 | 1.42.1 |
| tximport | 1.30.0 |
| clusterProfiler | 4.10.1 |
| R | 4.3.0 |

Reference genome: *Danio rerio* GRCz11 (GCF_000002035.6)

---

## Reference

Xue, Y.-H., et al. (2021). The time-dependent variations of zebrafish intestine and gill after polyethylene microplastics exposure. *Ecotoxicology*, 30(10), 1997–2010. https://doi.org/10.1007/s10646-021-02469-4
