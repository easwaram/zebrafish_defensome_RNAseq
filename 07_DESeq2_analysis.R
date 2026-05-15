# 07_DESeq2_analysis.R
# Differential gene expression analysis: gut vs gill in Danio rerio
# Input: Salmon quant files (tximport) and HTSeq count matrices
# Methods: tximport -> DESeq2, with apeglm LFC shrinkage
# Downstream: GO enrichment analysis of chemical defensome genes

# --- Install required packages (run once) ---
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2", "tximport", "tximeta",
  "GenomeInfoDb", "org.Dr.eg.db",
  "TxDb.Drerio.UCSC.danRer11.refGene",
  "clusterProfiler", "AnnotationDbi",
  "enrichplot", "apeglm"
))

install.packages(c("readr", "RColorBrewer", "gplots", "pheatmap", "ggupset", "tidyverse"))

# --- Libraries ---
library(DESeq2)
library(tximport)
library(TxDb.Drerio.UCSC.danRer11.refGene)
library(readr)
library(RColorBrewer)
library(gplots)
library(pheatmap)
library(tidyverse)
library(clusterProfiler)
library(org.Dr.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(ggupset)
library(apeglm)

# --- tximport: import Salmon quantification files ---

setwd("/path/to/GRCz11_salmon_quant")  # UPDATE this path

samples <- c(
  "SRR12614064_gut_C3",
  "SRR12614073_gut_C2",
  "SRR12614074_gut_C1",
  "SRR12876838_gill_C3",
  "SRR12876847_gill_C2",
  "SRR12876848_gill_C1"
)

quant_files <- file.path(dir(), "quant.sf")
names(quant_files) <- samples
stopifnot(all(file.exists(quant_files)))

# Build tx2gene mapping from TxDb
txdb <- TxDb.Drerio.UCSC.danRer11.refGene
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, keys = k, columns = "GENEID", keytype = "TXNAME")

# Import Salmon counts via tximport
txi <- tximport(
  quant_files,
  type = "salmon",
  tx2gene = tx2gene,
  ignoreTxVersion = TRUE  # removes version suffix from transcript IDs to match tx2gene
)

# Quick QC
cor_vals <- cor(txi$counts[, 1:6])
hist(cor_vals, main = "Sample Correlation Distribution")
pairs(log2(txi$counts[, 1:6]), pch = 20, lower.panel = NULL, col = rgb(0, 0, 1, 0.5))

# --- DESeq2 on Salmon tximport counts ---

tissue_design <- data.frame(
  row.names = samples,
  background = c(rep("gut", 3), rep("gill", 3))
)

dds_txi <- DESeqDataSetFromTximport(txi, colData = tissue_design, design = ~background)

# PCA
rlog_txi <- rlog(dds_txi, blind = TRUE)
plotPCA(rlog_txi, intgroup = "background", ntop = 2000)

# Run DESeq2
dds_txi2 <- DESeq(dds_txi)
plotDispEsts(dds_txi2, legend = FALSE)

# LFC shrinkage with apeglm
res_txi <- lfcShrink(dds_txi2, coef = "background_gut_vs_gill", type = "apeglm")
summary(res_txi)
plotMA(res_txi, ylim = c(-5, 5))

# --- filter for chemical defensome genes ---

defense_genes <- readLines("/path/to/ZFIN_genes/Complete lists/complete_ncbi_gene_ids")  # UPDATE

# Subset DESeq2 results to defensome genes
filtered_txi <- res_txi[defense_genes, ]
filtered_txi <- na.omit(filtered_txi)

# Re-adjust p-values for the defensome gene subset
filtered_txi$padj <- p.adjust(filtered_txi$pvalue, method = "BH")
summary(filtered_txi)

# Top DEGs
top_genes <- head(filtered_txi[order(filtered_txi$padj), ], 10)
print(top_genes)
write.csv(top_genes, file = "txi_top_genes.csv")

# --- visualizations ---

# Heatmap of top expressed genes
select_genes <- order(rowMeans(counts(dds_txi2, normalized = TRUE)), decreasing = TRUE)[1:6000]
df_annot <- as.data.frame(colData(dds_txi2)["background"])
ntd <- normTransform(dds_txi2)

pheatmap(
  assay(ntd)[select_genes, ],
  cluster_rows = FALSE,
  show_rownames = FALSE,
  cluster_cols = FALSE,
  annotation_col = df_annot
)

# Per-gene count plots (examples: abcg2c, cyp3c1)
plotCounts(dds_txi2, gene = "569858", intgroup = "background", pch = 6, col = "red")  # abcg2c
plotCounts(dds_txi2, gene = "324340", intgroup = "background", pch = 6, col = "red")  # cyp3c1

# --- GO enrichment analysis ---

sig_genes_txi <- filtered_txi[filtered_txi$padj < 0.05 & filtered_txi$baseMean > 50, ]
genes_up <- rownames(sig_genes_txi[sig_genes_txi$log2FoldChange > 0.5, ])

GO_results <- enrichGO(
  gene = genes_up,
  OrgDb = "org.Dr.eg.db",
  keyType = "ENTREZID",
  ont = "BP"
)

# Bar plot of GO terms
fit <- barplot(GO_results, showCategory = 20)
png("GO_barplot_txi_up.png", res = 250, width = 1000, height = 1200)
print(fit)
dev.off()

# Directional gene lists for Venn diagram
txi_up <- filtered_txi[filtered_txi$log2FoldChange > 0 & filtered_txi$padj < 0.1, ]
txi_down <- filtered_txi[filtered_txi$log2FoldChange < 0 & filtered_txi$padj < 0.1, ]
write.csv(rownames(txi_up), file = "txi_up_genes.csv")
write.csv(rownames(txi_down), file = "txi_down_genes.csv")

# --- HTSeq count matrix (alternative to Salmon) ---

setwd("/path/to/htseq_counts")  # UPDATE

read_sample <- function(sample_name) {
  fname <- paste0(sample_name, "_htseq_counts.txt")
  read.delim(fname, col.names = c("gene", "count"), sep = "\t",
             colClasses = c("character", "numeric"))
}

all_data <- read_sample(samples[1])
for (s in samples[-1]) {
  temp <- read_sample(s)
  all_data <- cbind(all_data, temp$count)
}
colnames(all_data)[2:ncol(all_data)] <- samples

# Remove HTSeq summary rows at bottom
all_data <- all_data[1:(nrow(all_data) - 5), ]

raw_counts <- all_data[, 2:ncol(all_data)]
rownames(raw_counts) <- all_data$gene

# Build DESeq2 object from HTSeq matrix
dds_mat <- DESeqDataSetFromMatrix(
  countData = raw_counts,
  colData = tissue_design,
  design = ~background
)

rlog_mat <- rlog(dds_mat, blind = TRUE)
plotPCA(rlog_mat, intgroup = "background", ntop = 2000)

dds_mat2 <- DESeq(dds_mat)
res_mat <- lfcShrink(dds_mat2, coef = "background_gut_vs_gill", type = "apeglm")
summary(res_mat)
plotMA(res_mat, ylim = c(-5, 5))

# Filter for defensome genes and re-adjust p-values
filtered_mat <- res_mat[defense_genes, ]
filtered_mat <- na.omit(filtered_mat)
filtered_mat$padj <- p.adjust(filtered_mat$pvalue, method = "BH")

top_genes_mat <- head(filtered_mat[order(filtered_mat$padj), ], 10)
write.csv(top_genes_mat, file = "mat_top_genes.csv")
