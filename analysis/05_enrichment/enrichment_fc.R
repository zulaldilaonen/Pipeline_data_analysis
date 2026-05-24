#!/usr/bin/env Rscript

# Purpose:
#   Perform ORA and GSEA enrichment analysis for frontal cortex DESeq2 results.
#
# How to run:
#   cd ~/Pipeline_data_analysis/analysis
#   Rscript 05_enrichment/enrichment_fc.R

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
})

# Input/output paths
input_file <- "03_deseq2/results_diagnosis_in_FC.csv"
out_dir <- "05_enrichment"

# Read DESeq2 results
res <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)

# The first column is usually the gene ID column written from rownames
colnames(res)[1] <- "gene_id"

# Remove version numbers from Ensembl IDs: ENSG000001234.5 -> ENSG000001234
res$ENSEMBL <- sub("\\..*$", "", res$gene_id)

# Remove rows with missing p-values/statistics
res <- res[!is.na(res$pvalue) & !is.na(res$stat), ]

# Remove duplicated Ensembl IDs if any
res <- res[!duplicated(res$ENSEMBL), ]

# -----------------------------
# 1) ORA using significant genes
# -----------------------------
sig_res <- res[!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 1, ]

# Convert significant genes to Entrez IDs
sig_map <- bitr(sig_res$ENSEMBL,
                fromType = "ENSEMBL",
                toType = c("ENTREZID", "SYMBOL"),
                OrgDb = org.Hs.eg.db)

sig_entrez <- unique(sig_map$ENTREZID)

# Background universe (recommended)
bg_map <- bitr(res$ENSEMBL,
               fromType = "ENSEMBL",
               toType = c("ENTREZID"),
               OrgDb = org.Hs.eg.db)

bg_entrez <- unique(bg_map$ENTREZID)

# ORA GO
ora_go <- NULL
if (length(sig_entrez) > 0) {
  ora_go <- enrichGO(
    gene = sig_entrez,
    universe = bg_entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
}

# ORA KEGG
ora_kegg <- NULL
if (length(sig_entrez) > 0) {
  ora_kegg <- enrichKEGG(
    gene = sig_entrez,
    universe = bg_entrez,
    organism = "hsa",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )
}

# Save ORA tables
if (!is.null(ora_go) && nrow(as.data.frame(ora_go)) > 0) {
  write.csv(as.data.frame(ora_go), file.path(out_dir, "ora_go_fc.csv"), row.names = FALSE)
  png(file.path(out_dir, "ora_go_fc_dotplot.png"), width = 1800, height = 1400, res = 200)
  print(dotplot(ora_go, showCategory = 20) + ggtitle("ORA GO BP - FC"))
  dev.off()
}

if (!is.null(ora_kegg) && nrow(as.data.frame(ora_kegg)) > 0) {
  write.csv(as.data.frame(ora_kegg), file.path(out_dir, "ora_kegg_fc.csv"), row.names = FALSE)
  png(file.path(out_dir, "ora_kegg_fc_dotplot.png"), width = 1800, height = 1400, res = 200)
  print(dotplot(ora_kegg, showCategory = 20) + ggtitle("ORA KEGG - FC"))
  dev.off()
}

# -----------------------------
# 2) GSEA using all ranked genes
# -----------------------------
gsea_map <- bitr(res$ENSEMBL,
                 fromType = "ENSEMBL",
                 toType = c("ENTREZID"),
                 OrgDb = org.Hs.eg.db)

# Keep only mapped genes
res_gsea <- res[res$ENSEMBL %in% gsea_map$ENSEMBL, ]

# Merge to get ENTREZID
res_gsea <- merge(res_gsea, gsea_map, by = "ENSEMBL")

# If one ENTREZ maps multiple times, keep the row with highest absolute stat
res_gsea <- res_gsea[order(abs(res_gsea$stat), decreasing = TRUE), ]
res_gsea <- res_gsea[!duplicated(res_gsea$ENTREZID), ]

# Ranked gene list
gene_list <- res_gsea$stat
names(gene_list) <- res_gsea$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

# GSEA GO
gsea_go <- gseGO(
  geneList = gene_list,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE
)

# GSEA KEGG
gsea_kegg <- gseKEGG(
  geneList = gene_list,
  organism = "hsa",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE
)

# Save GSEA tables
if (!is.null(gsea_go) && nrow(as.data.frame(gsea_go)) > 0) {
  write.csv(as.data.frame(gsea_go), file.path(out_dir, "gsea_go_fc.csv"), row.names = FALSE)
  png(file.path(out_dir, "gsea_go_fc_dotplot.png"), width = 1800, height = 1400, res = 200)
  print(dotplot(gsea_go, showCategory = 20) + ggtitle("GSEA GO BP - FC"))
  dev.off()
}

if (!is.null(gsea_kegg) && nrow(as.data.frame(gsea_kegg)) > 0) {
  write.csv(as.data.frame(gsea_kegg), file.path(out_dir, "gsea_kegg_fc.csv"), row.names = FALSE)
  png(file.path(out_dir, "gsea_kegg_fc_dotplot.png"), width = 1800, height = 1400, res = 200)
  print(dotplot(gsea_kegg, showCategory = 20) + ggtitle("GSEA KEGG - FC"))
  dev.off()
}

# Summary file
sink(file.path(out_dir, "enrichment_summary.txt"))
cat("Enrichment Summary - FC\n")
cat("=======================\n\n")
cat("Input file:", input_file, "\n")
cat("Total genes after cleanup:", nrow(res), "\n")
cat("Significant genes for ORA (padj < 0.05 and |log2FC| > 1):", nrow(sig_res), "\n")
cat("Mapped significant Entrez IDs:", length(sig_entrez), "\n")
cat("Mapped background Entrez IDs:", length(bg_entrez), "\n")
cat("Genes in GSEA ranked list:", length(gene_list), "\n\n")

if (!is.null(ora_go)) cat("ORA GO terms:", nrow(as.data.frame(ora_go)), "\n")
if (!is.null(ora_kegg)) cat("ORA KEGG terms:", nrow(as.data.frame(ora_kegg)), "\n")
if (!is.null(gsea_go)) cat("GSEA GO terms:", nrow(as.data.frame(gsea_go)), "\n")
if (!is.null(gsea_kegg)) cat("GSEA KEGG terms:", nrow(as.data.frame(gsea_kegg)), "\n")
sink()

cat("Enrichment analysis completed successfully.\n")
