library(DESeq2)
library(ggplot2)
library(pheatmap)

# Load data
dds <- readRDS("03_deseq2/dds.rds")
meta <- read.csv("00_metadata/metadata_filtered.csv", row.names=1)

# PCA
vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("diagnosis","brain_region"), returnData=TRUE)

p <- ggplot(pcaData, aes(PC1, PC2, color=diagnosis, shape=brain_region)) +
  geom_point(size=3) +
  theme_minimal()

ggsave("04_visualization/pca_plot.png", p)

# Volcano function
plot_volcano <- function(file, outname) {
  df <- read.csv(file)
  
  df$significant <- "no"
  df$significant[df$padj < 0.05 & abs(df$log2FoldChange) > 1] <- "yes"
  
  p <- ggplot(df, aes(x=log2FoldChange, y=-log10(padj), color=significant)) +
    geom_point(alpha=0.6) +
    theme_minimal() +
    scale_color_manual(values=c("grey","red"))
  
  ggsave(outname, p)
}

# Volcano plots
plot_volcano("03_deseq2/results_diagnosis_in_FC.csv", "04_visualization/volcano_fc.png")
plot_volcano("03_deseq2/results_diagnosis_in_CB.csv", "04_visualization/volcano_cb.png")

# Heatmap (top 50 genes FC)
res_fc <- read.csv("03_deseq2/results_diagnosis_in_FC.csv", row.names=1)

res_fc <- res_fc[order(res_fc$padj), ]
top_genes <- rownames(res_fc)[1:50]

mat <- assay(vsd)[top_genes, ]

# scale rows
mat <- t(scale(t(mat)))

pheatmap(mat,
         annotation_col = meta[,c("diagnosis","brain_region")],
         filename="04_visualization/heatmap_top_genes.png")

print("Visualization completed")
