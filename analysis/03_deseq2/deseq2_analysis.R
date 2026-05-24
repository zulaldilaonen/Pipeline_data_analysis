library(DESeq2)

# Load data
counts <- read.table("01_counts/merged_counts_filtered.tsv", header=TRUE, row.names=1, check.names=FALSE)
meta <- read.csv("00_metadata/metadata_filtered.csv", row.names=1, check.names=FALSE)

# Match sample order
meta <- meta[colnames(counts), ]

# Convert variables
meta$sex <- factor(meta$sex)
meta$brain_region <- factor(meta$brain_region, levels = c("CB", "FC"))
meta$diagnosis <- factor(meta$diagnosis, levels = c("Control", "AD"))

# Center numeric covariates
meta$age_centered <- scale(meta$age, center = TRUE, scale = FALSE)
meta$pmi_centered <- scale(meta$pmi, center = TRUE, scale = FALSE)

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta,
  design = ~ age_centered + sex + pmi_centered + brain_region + diagnosis + brain_region:diagnosis
)

# Optional low-count filtering
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Run DESeq2
dds <- DESeq(dds)

# Save coefficient names
coef_names <- resultsNames(dds)
write.table(coef_names,
            file = "03_deseq2/results_names.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

# Main diagnosis effect at reference brain region (CB)
res_diagnosis_cb <- results(dds, name = "diagnosis_AD_vs_Control")
write.csv(as.data.frame(res_diagnosis_cb), "03_deseq2/results_diagnosis_in_CB.csv")

# Interaction term: whether AD effect differs in FC vs CB
interaction_name <- grep("brain_regionFC.*diagnosisAD|diagnosisAD.*brain_regionFC",
                         coef_names, value = TRUE)

if (length(interaction_name) == 1) {
  res_interaction <- results(dds, name = interaction_name)
  write.csv(as.data.frame(res_interaction), "03_deseq2/results_interaction_FC_vs_CB.csv")
}

# AD effect in FC = diagnosis main effect + interaction
if (length(interaction_name) == 1) {
  res_diagnosis_fc <- results(
    dds,
    contrast = list(c("diagnosis_AD_vs_Control", interaction_name))
  )
  write.csv(as.data.frame(res_diagnosis_fc), "03_deseq2/results_diagnosis_in_FC.csv")
}

# Shrink only if the coefficient exists
if ("diagnosis_AD_vs_Control" %in% coef_names) {
  res_shrink_cb <- lfcShrink(dds, coef = "diagnosis_AD_vs_Control", type = "apeglm")
  write.csv(as.data.frame(res_shrink_cb), "03_deseq2/results_diagnosis_in_CB_shrunken.csv")
}

# Save DDS object
saveRDS(dds, file = "03_deseq2/dds.rds")

print("DESeq2 analysis completed successfully.")
print("Available coefficient names:")
print(coef_names)
