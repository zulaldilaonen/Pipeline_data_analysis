#!/usr/bin/env Rscript
# =============================================================================
# figure_generation_script.R
# -----------------------------------------------------------------------------
# Purpose
#   Tek bir komut ile çalıştırıldığında, mevcut DESeq2 / enrichment çıktılarını
#   bozmadan, akademik tez formatına uygun publication-ready görseller üretir.
#
#   Bu script SADECE görselleri ve not dosyalarını yeniden üretir. Raw count
#   matrix, metadata, differential expression CSV ve enrichment CSV
#   dosyalarını DEĞİŞTİRMEZ.
#
# How to run
#   cd ~/Pipeline_data_analysis/analysis
#   Rscript figure_generation_script.R
#
# Beklenen klasör yapısı (orijinal pipeline ile aynı):
#   00_metadata/metadata_filtered.csv
#   01_counts/merged_counts_filtered.tsv     (PCA/heatmap için gerekli olabilir)
#   03_deseq2/dds.rds                        (PCA/heatmap için gerekli)
#   03_deseq2/results_diagnosis_in_FC.csv
#   03_deseq2/results_diagnosis_in_CB.csv
#   05_enrichment/ora_go_fc.csv
#   05_enrichment/gsea_go_fc.csv
#   05_enrichment/gsea_kegg_fc.csv
#
# Çıktılar (yeni dosyalar, eski dosyalara dokunulmaz):
#   04_visualization/volcano_fc_publication.{png,pdf}
#   04_visualization/volcano_cb_publication.{png,pdf}
#   04_visualization/pca_plot_publication.{png,pdf}
#   04_visualization/heatmap_top_genes_publication.{png,pdf}
#   05_enrichment/ora_go_fc_dotplot_publication.{png,pdf}
#   05_enrichment/gsea_go_fc_dotplot_publication.{png,pdf}
#   05_enrichment/gsea_kegg_fc_dotplot_publication.{png,pdf}
#   04_visualization/deseq2_version.txt
#   04_visualization/methods_notes.txt
#   04_visualization/figure_notes.txt
# =============================================================================


# -----------------------------------------------------------------------------
# 0. KONFİGÜRASYON
# -----------------------------------------------------------------------------
# Yollar
METADATA_FILE      <- "00_metadata/metadata_filtered.csv"
DDS_FILE           <- "03_deseq2/dds.rds"
RES_FC_FILE        <- "03_deseq2/results_diagnosis_in_FC.csv"
RES_CB_FILE        <- "03_deseq2/results_diagnosis_in_CB.csv"

ORA_GO_FC_FILE     <- "05_enrichment/ora_go_fc.csv"
GSEA_GO_FC_FILE    <- "05_enrichment/gsea_go_fc.csv"
GSEA_KEGG_FC_FILE  <- "05_enrichment/gsea_kegg_fc.csv"

VIS_DIR    <- "04_visualization"
ENRICH_DIR <- "05_enrichment"

# Significance thresholds (tek noktada)
PADJ_CUTOFF   <- 0.05
LFC_CUTOFF    <- 1
N_LABEL_GENES <- 12          # volcano plot'ta etiketlenecek gen sayısı
N_HEATMAP     <- 50          # heatmap'te gösterilecek top gen sayısı
N_DOTPLOT     <- 15          # dotplot'ta gösterilecek top pathway sayısı

# Üretilen / atlanan görsel kayıtları
.produced <- character(0)
.skipped  <- character(0)
.notes    <- character(0)
add_note  <- function(msg) .notes <<- c(.notes, msg)
mark_done <- function(...) .produced <<- c(.produced, c(...))
mark_skip <- function(name, reason) .skipped <<- c(.skipped, paste0(name, "  --  ", reason))

# Klasörleri oluştur
dir.create(VIS_DIR,    showWarnings = FALSE, recursive = TRUE)
dir.create(ENRICH_DIR, showWarnings = FALSE, recursive = TRUE)


# -----------------------------------------------------------------------------
# 1. PAKET KONTROLÜ
# -----------------------------------------------------------------------------
# Zorunlu paketler (yoksa açık hata)
required_pkgs <- c("DESeq2", "ggplot2", "ggrepel", "pheatmap",
                   "dplyr", "stringr")

missing_pkgs <- required_pkgs[!vapply(required_pkgs,
                                      requireNamespace, logical(1),
                                      quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Eksik paket(ler): ", paste(missing_pkgs, collapse = ", "), "\n",
    "Kurulum:\n",
    "  CRAN:        install.packages(c('ggplot2','ggrepel','pheatmap',",
    "'dplyr','readr','stringr','tibble'))\n",
    "  Bioconductor: BiocManager::install(c('DESeq2'))\n"
  )
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(dplyr)
  library(stringr)
})

# Opsiyonel: ENSEMBL -> SYMBOL haritası için
has_orgdb <- requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
             requireNamespace("AnnotationDbi", quietly = TRUE)


# -----------------------------------------------------------------------------
# 2. ORTAK AKADEMİK TEMA
# -----------------------------------------------------------------------------
theme_pub <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major   = element_line(colour = "grey92", linewidth = 0.3),
      plot.title         = element_text(face = "bold", hjust = 0.5,
                                        size = base_size + 2),
      plot.subtitle      = element_text(hjust = 0.5, size = base_size),
      axis.title         = element_text(face = "bold"),
      axis.text          = element_text(colour = "black"),
      legend.title       = element_text(face = "bold"),
      legend.position    = "right",
      legend.background  = element_blank(),
      strip.background   = element_rect(fill = "grey95", colour = NA),
      strip.text         = element_text(face = "bold")
    )
}

# Tutarlı renkler
COL_UP    <- "#C0392B"   # kırmızı  - upregulated
COL_DOWN  <- "#1F618D"   # mavi     - downregulated
COL_NS    <- "grey70"    # gri      - not significant
COL_GROUP <- c("Control" = "#2E86AB", "AD" = "#E63946")


# -----------------------------------------------------------------------------
# 3. YARDIMCI FONKSİYONLAR
# -----------------------------------------------------------------------------

# 300 dpi PNG + PDF kaydet
save_pub <- function(plot_obj, base_path, width = 8, height = 6) {
  ggsave(paste0(base_path, ".png"), plot_obj,
         width = width, height = height, dpi = 300, units = "in")
  ggsave(paste0(base_path, ".pdf"), plot_obj,
         width = width, height = height, units = "in")
  mark_done(paste0(base_path, ".png"), paste0(base_path, ".pdf"))
}

# Volcano için zorunlu kolonları algıla / standardize et
prepare_de_table <- function(df, source_label) {
  # Birinci kolon genelde rowname olarak yazılmış (boş header)
  if (colnames(df)[1] %in% c("", "X", "X1", "gene", "gene_id", "ID")) {
    colnames(df)[1] <- "gene_id"
  } else if (!"gene_id" %in% colnames(df)) {
    df$gene_id <- rownames(df)
  }
  required <- c("gene_id", "log2FoldChange", "padj")
  miss <- setdiff(required, colnames(df))
  if (length(miss) > 0) {
    stop(sprintf("[%s] Eksik kolon(lar): %s",
                 source_label, paste(miss, collapse = ", ")))
  }
  df
}

# ENSEMBL'i sürümsüz ve sembolle birlikte getir
add_gene_symbols <- function(df) {
  df$ensembl_clean <- sub("\\..*$", "", df$gene_id)
  if (has_orgdb) {
    sym <- tryCatch(
      AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                            keys      = df$ensembl_clean,
                            column    = "SYMBOL",
                            keytype   = "ENSEMBL",
                            multiVals = "first"),
      error = function(e) NULL
    )
    df$symbol <- if (!is.null(sym)) unname(sym) else NA_character_
  } else {
    df$symbol <- NA_character_
  }
  # Etiketlerde kullanılacak görüntü ismi
  df$label <- ifelse(is.na(df$symbol) | df$symbol == "",
                     df$ensembl_clean, df$symbol)
  df
}


# -----------------------------------------------------------------------------
# 4. METHODS / VERSION NOTLARI
# -----------------------------------------------------------------------------
deseq2_version <- as.character(packageVersion("DESeq2"))
writeLines(deseq2_version, file.path(VIS_DIR, "deseq2_version.txt"))
mark_done(file.path(VIS_DIR, "deseq2_version.txt"))

methods_lines <- c(
  "Methods notes — DESeq2 differential expression and visualization",
  paste(rep("=", 70), collapse = ""),
  "",
  paste0("DESeq2 version: ", deseq2_version),
  paste0("R version:      ", paste(R.version$major, R.version$minor, sep = ".")),
  paste0("Run date:       ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "Normalization",
  "  DESeq2 uses the median-of-ratios size factor normalization",
  "  (estimateSizeFactors). Library size and composition biases are",
  "  corrected by computing per-sample size factors from the geometric",
  "  mean of gene counts across samples.",
  "",
  "Dispersion and testing",
  "  Gene-wise dispersions are estimated with the Cox-Reid adjusted",
  "  maximum likelihood, shrunk toward a fitted mean-dispersion trend.",
  "  Differential expression is tested with the Wald test on the",
  "  shrunken dispersion estimates (default DESeq2 workflow).",
  "",
  "Multiple-testing correction",
  "  Adjusted p-values reported in results() use the Benjamini-Hochberg",
  "  procedure (pAdjustMethod = 'BH'), which is the DESeq2 default.",
  "",
  "Variance-stabilizing transformation",
  "  PCA and heatmap use vst(dds, blind = FALSE) to obtain homoscedastic,",
  "  log-scale expression values while accounting for the experimental",
  "  design.",
  "",
  "Significance thresholds (applied uniformly in this script)",
  paste0("  Adjusted p-value (padj) cutoff: ", PADJ_CUTOFF),
  paste0("  log2 fold-change cutoff:        |log2FC| > ", LFC_CUTOFF),
  "  Gene classification:",
  paste0("    Upregulated:   padj < ", PADJ_CUTOFF,
         " AND log2FoldChange >  ", LFC_CUTOFF),
  paste0("    Downregulated: padj < ", PADJ_CUTOFF,
         " AND log2FoldChange < -", LFC_CUTOFF),
  "    Not significant: all other genes",
  "",
  "Heatmap gene selection",
  paste0("  Top ", N_HEATMAP,
         " genes ranked by ascending padj in the FC contrast,"),
  "  visualized on vst-transformed counts with row-wise z-score scaling.",
  "",
  "Sample ordering in heatmap",
  "  Samples are ordered so that Control samples are grouped on one side",
  "  and AD samples on the other; within each diagnosis group, samples",
  "  are ordered by brain region.",
  "",
  "Enrichment analysis (re-using existing CSV outputs)",
  "  ORA/GSEA computed previously with clusterProfiler, p.adjust = 'BH'."
)
writeLines(methods_lines, file.path(VIS_DIR, "methods_notes.txt"))
mark_done(file.path(VIS_DIR, "methods_notes.txt"))


# -----------------------------------------------------------------------------
# 5. VOLCANO PLOT
# -----------------------------------------------------------------------------
make_volcano <- function(res_path, out_base, title_text) {
  if (!file.exists(res_path)) {
    mark_skip(out_base, sprintf("Girdi dosyası bulunamadı: %s", res_path))
    return(invisible(NULL))
  }

  df <- read.csv(res_path, stringsAsFactors = FALSE, check.names = FALSE)
  df <- prepare_de_table(df, basename(res_path))
  df <- add_gene_symbols(df)

  # NA padj olan satırları (DESeq2 independent filtering kaynaklı) at
  df <- df[!is.na(df$padj) & !is.na(df$log2FoldChange), ]

  # Sınıflandırma
  df$category <- "Not significant"
  df$category[df$padj < PADJ_CUTOFF & df$log2FoldChange >  LFC_CUTOFF] <- "Upregulated"
  df$category[df$padj < PADJ_CUTOFF & df$log2FoldChange < -LFC_CUTOFF] <- "Downregulated"
  df$category <- factor(df$category,
                        levels = c("Upregulated", "Downregulated", "Not significant"))

  # En anlamlı genlerden etiket için ilk N tane (yalnızca up/down arasından)
  sig_df <- df[df$category != "Not significant", ]
  sig_df <- sig_df[order(sig_df$padj), ]
  top_df <- head(sig_df, N_LABEL_GENES)

  # Sonsuz -log10(padj) değerlerini sınırla (görsellik için)
  max_neglog <- max(-log10(df$padj[df$padj > 0]), na.rm = TRUE)
  cap_val <- max_neglog * 1.05
  df$neglog10padj     <- pmin(-log10(df$padj), cap_val)
  top_df$neglog10padj <- pmin(-log10(top_df$padj), cap_val)

  p <- ggplot(df, aes(x = log2FoldChange, y = neglog10padj, colour = category)) +
    geom_point(alpha = 0.55, size = 1.4) +
    geom_hline(yintercept = -log10(PADJ_CUTOFF),
               linetype = "dashed", colour = "grey30", linewidth = 0.4) +
    geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF),
               linetype = "dashed", colour = "grey30", linewidth = 0.4) +
    scale_colour_manual(values = c("Upregulated"     = COL_UP,
                                   "Downregulated"   = COL_DOWN,
                                   "Not significant" = COL_NS),
                        name = "Regulation") +
    labs(
      title    = title_text,
      subtitle = sprintf("padj < %.2g  |  |log2FC| > %g",
                         PADJ_CUTOFF, LFC_CUTOFF),
      x = expression(log[2]~"fold change"),
      y = expression(-log[10]~"(adjusted p-value)")
    ) +
    theme_pub()

  if (nrow(top_df) > 0) {
    p <- p +
      geom_text_repel(
        data            = top_df,
        aes(label = label),
        size            = 3.2,
        colour          = "black",
        box.padding     = 0.5,
        point.padding   = 0.3,
        segment.colour  = "grey40",
        segment.size    = 0.3,
        max.overlaps    = Inf,
        min.segment.length = 0,
        show.legend     = FALSE
      )
  }

  save_pub(p, file.path(VIS_DIR, out_base), width = 8, height = 6.5)

  add_note(sprintf(
    "Volcano (%s): %d up, %d down, %d not significant; %d genes labelled.",
    out_base,
    sum(df$category == "Upregulated"),
    sum(df$category == "Downregulated"),
    sum(df$category == "Not significant"),
    nrow(top_df)
  ))
}

make_volcano(RES_FC_FILE, "volcano_fc_publication",
             "Volcano plot — AD vs. Control in Frontal Cortex (FC)")
make_volcano(RES_CB_FILE, "volcano_cb_publication",
             "Volcano plot — AD vs. Control in Cerebellum (CB)")


# -----------------------------------------------------------------------------
# 6. PCA PLOT
# -----------------------------------------------------------------------------
dds <- NULL
vsd <- NULL
if (file.exists(DDS_FILE)) {
  dds <- readRDS(DDS_FILE)
  vsd <- vst(dds, blind = FALSE)
} else {
  add_note("PCA ve heatmap için dds.rds bulunamadı — bu görseller atlandı.")
}

if (!is.null(vsd)) {
  intgroup <- intersect(c("diagnosis", "brain_region"), colnames(colData(vsd)))
  if (length(intgroup) == 0) {
    mark_skip("pca_plot_publication",
              "colData içinde diagnosis/brain_region bulunamadı.")
  } else {
    pca_data <- plotPCA(vsd, intgroup = intgroup, returnData = TRUE)
    pct_var  <- round(100 * attr(pca_data, "percentVar"), 1)

    pca_data$diagnosis <- factor(pca_data$diagnosis,
                                 levels = c("Control", "AD"))

    p_pca <- ggplot(pca_data,
                    aes(x = PC1, y = PC2,
                        colour = diagnosis,
                        shape  = brain_region)) +
      geom_point(size = 3, alpha = 0.85, stroke = 0.6) +
      scale_colour_manual(values = COL_GROUP, name = "Diagnosis") +
      scale_shape_manual(values = c("CB" = 16, "FC" = 17),
                         name   = "Brain region") +
      labs(
        title = "PCA Plot of Variance-Stabilized Expression Data",
        x     = sprintf("PC1 (%.1f%% variance)", pct_var[1]),
        y     = sprintf("PC2 (%.1f%% variance)", pct_var[2])
      ) +
      theme_pub() +
      guides(colour = guide_legend(order = 1),
             shape  = guide_legend(order = 2))

    save_pub(p_pca, file.path(VIS_DIR, "pca_plot_publication"),
             width = 8, height = 6)
    add_note(sprintf(
      "PCA: vst(dds, blind = FALSE); %d sample, PC1 = %.1f%%, PC2 = %.1f%%.",
      nrow(pca_data), pct_var[1], pct_var[2]
    ))
  }
} else {
  mark_skip("pca_plot_publication", "dds.rds yok.")
}


# -----------------------------------------------------------------------------
# 7. HEATMAP (top genes, FC kontrastı, gruplara göre sıralı)
# -----------------------------------------------------------------------------
if (!is.null(vsd) && file.exists(RES_FC_FILE) && file.exists(METADATA_FILE)) {

  meta_full <- read.csv(METADATA_FILE, stringsAsFactors = FALSE,
                        check.names = FALSE)
  if (!"sample" %in% colnames(meta_full)) {
    if (colnames(meta_full)[1] == "") {
      colnames(meta_full)[1] <- "sample"
    } else {
      meta_full$sample <- rownames(meta_full)
    }
  }
  rownames(meta_full) <- meta_full$sample

  res_fc <- read.csv(RES_FC_FILE, stringsAsFactors = FALSE,
                     check.names = FALSE)
  res_fc <- prepare_de_table(res_fc, basename(RES_FC_FILE))

  # En anlamlı padj'a göre top N gen
  res_fc <- res_fc[!is.na(res_fc$padj), ]
  res_fc <- res_fc[order(res_fc$padj), ]
  top_genes <- res_fc$gene_id[seq_len(min(N_HEATMAP, nrow(res_fc)))]

  # vsd üzerinden sadece gerçekten var olanları çek
  vst_mat   <- assay(vsd)
  top_genes <- intersect(top_genes, rownames(vst_mat))

  if (length(top_genes) >= 2) {
    mat <- vst_mat[top_genes, , drop = FALSE]
    mat <- t(scale(t(mat)))                  # row-wise z-score

    # Satır isimlerini SYMBOL'a çevir (varsa)
    if (has_orgdb) {
      ens_clean <- sub("\\..*$", "", rownames(mat))
      sym <- tryCatch(
        AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                              keys      = ens_clean,
                              column    = "SYMBOL",
                              keytype   = "ENSEMBL",
                              multiVals = "first"),
        error = function(e) NULL
      )
      if (!is.null(sym)) {
        new_names <- ifelse(is.na(sym) | sym == "", ens_clean, unname(sym))
        # Eşsiz yap
        new_names <- make.unique(new_names)
        rownames(mat) <- new_names
      }
    }

    # Annotation tablosu
    ann_cols <- intersect(c("diagnosis", "brain_region"), colnames(meta_full))
    ann_df   <- meta_full[colnames(mat), ann_cols, drop = FALSE]
    ann_df$diagnosis <- factor(ann_df$diagnosis,
                               levels = c("Control", "AD"))
    if ("brain_region" %in% colnames(ann_df)) {
      ann_df$brain_region <- factor(ann_df$brain_region,
                                    levels = c("CB", "FC"))
    }

    # Örnek sırası: Control önce, AD sonra; içinde brain_region'a göre alt-sıralama
    order_idx <- order(ann_df$diagnosis,
                       if ("brain_region" %in% colnames(ann_df))
                         ann_df$brain_region else rep(1, nrow(ann_df)),
                       colnames(mat))
    mat    <- mat[, order_idx, drop = FALSE]
    ann_df <- ann_df[order_idx, , drop = FALSE]

    ann_colors <- list(
      diagnosis    = COL_GROUP,
      brain_region = c("CB" = "#88B04B", "FC" = "#F4A261")
    )

    # Sütun isimleri çoğu zaman çok kalabalık — gizle, annotation grupları gösterir
    show_colnames <- (ncol(mat) <= 30)

    # Heatmap'i grafik device üzerinden çizmek için pheatmap'i bir gridob olarak alıp
    # hem PNG hem PDF üretiyoruz.
    draw_hm <- function() {
      pheatmap(
        mat,
        annotation_col    = ann_df,
        annotation_colors = ann_colors,
        cluster_cols      = FALSE,                 # grup sırasını korusun
        cluster_rows      = TRUE,
        show_colnames     = show_colnames,
        show_rownames     = TRUE,
        fontsize          = 8,
        fontsize_row      = 7,
        fontsize_col      = 6,
        treeheight_row    = 20,
        treeheight_col    = 0,
        border_color      = NA,
        color             = colorRampPalette(
                              c("#1F618D", "white", "#C0392B"))(100),
        main              = sprintf(
          "Top %d DE genes (FC contrast) — vst, row z-score",
          length(top_genes)),
        silent            = TRUE
      )
    }

    hm_obj <- draw_hm()

    # PNG
    png(file.path(VIS_DIR, "heatmap_top_genes_publication.png"),
        width = 10, height = 11, units = "in", res = 300)
    grid::grid.newpage(); grid::grid.draw(hm_obj$gtable)
    dev.off()

    # PDF
    pdf(file.path(VIS_DIR, "heatmap_top_genes_publication.pdf"),
        width = 10, height = 11)
    grid::grid.newpage(); grid::grid.draw(hm_obj$gtable)
    dev.off()

    mark_done(
      file.path(VIS_DIR, "heatmap_top_genes_publication.png"),
      file.path(VIS_DIR, "heatmap_top_genes_publication.pdf")
    )
    add_note(sprintf(
      "Heatmap: top %d genes by FC padj; column clustering disabled; samples ordered Control → AD (by brain_region within).",
      length(top_genes)
    ))
  } else {
    mark_skip("heatmap_top_genes_publication",
              "vsd'de eşleşen top gene sayısı < 2.")
  }
} else {
  mark_skip("heatmap_top_genes_publication",
            "Gerekli dosyalar eksik (dds.rds / FC results / metadata).")
}


# -----------------------------------------------------------------------------
# 8. ENRICHMENT DOTPLOTLARI (ORA GO, GSEA GO, GSEA KEGG)
# -----------------------------------------------------------------------------

# clusterProfiler nesnesi olmadığı için dotplot()'u kullanamıyoruz;
# CSV'lerden ggplot ile yeniden çiziyoruz. Stil tüm dotplotlarda ortak.
plot_enrichment_dotplot <- function(csv_path, out_base, title_text,
                                    mode = c("ora", "gsea")) {
  mode <- match.arg(mode)
  if (!file.exists(csv_path)) {
    mark_skip(out_base, sprintf("Girdi dosyası bulunamadı: %s", csv_path))
    return(invisible(NULL))
  }
  df <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)

  needed <- c("Description", "p.adjust")
  miss <- setdiff(needed, colnames(df))
  if (length(miss) > 0) {
    mark_skip(out_base,
              sprintf("Eksik kolon(lar): %s", paste(miss, collapse = ", ")))
    return(invisible(NULL))
  }

  if (mode == "ora") {
    # GeneRatio "a/b" -> sayısal oran
    parse_ratio <- function(x) {
      vapply(strsplit(x, "/", fixed = TRUE), function(p) {
        if (length(p) == 2) as.numeric(p[1]) / as.numeric(p[2]) else NA_real_
      }, numeric(1))
    }
    if ("GeneRatio" %in% colnames(df)) {
      df$gene_ratio_num <- parse_ratio(df$GeneRatio)
    } else {
      df$gene_ratio_num <- NA_real_
    }
    df$x_axis_value <- df$gene_ratio_num
    df$size_value   <- if ("Count" %in% colnames(df)) df$Count else NA_real_
    x_label <- "Gene ratio"
    size_label <- "Count"
  } else {
    # GSEA
    df$x_axis_value <- if ("NES" %in% colnames(df)) df$NES else NA_real_
    df$size_value   <- if ("setSize" %in% colnames(df)) df$setSize else NA_real_
    x_label <- "Normalized enrichment score (NES)"
    size_label <- "Set size"
  }

  df <- df[!is.na(df$p.adjust), ]
  df <- df[order(df$p.adjust), ]
  df <- head(df, N_DOTPLOT)

  if (nrow(df) == 0) {
    mark_skip(out_base, "Çizilecek pathway kalmadı.")
    return(invisible(NULL))
  }

  # Uzun pathway isimlerini sar; sıralama p.adjust'a göre
  df$Description_wrap <- str_wrap(df$Description, width = 45)
  df$Description_wrap <- factor(df$Description_wrap,
                                levels = rev(df$Description_wrap))

  p <- ggplot(df, aes(x = x_axis_value, y = Description_wrap)) +
    geom_point(aes(size = size_value, colour = p.adjust)) +
    scale_colour_gradient(low = "#C0392B", high = "#1F618D",
                          name = "p.adjust",
                          guide = guide_colourbar(reverse = FALSE)) +
    scale_size_continuous(name = size_label, range = c(2, 7)) +
    labs(title = title_text, x = x_label, y = NULL) +
    theme_pub(base_size = 11) +
    theme(axis.text.y = element_text(size = 9, lineheight = 0.9))

  if (mode == "gsea") {
    p <- p + geom_vline(xintercept = 0, linetype = "dashed",
                        colour = "grey40", linewidth = 0.3)
  }

  save_pub(p, file.path(ENRICH_DIR, out_base),
           width = 9, height = 6 + 0.18 * nrow(df))
  add_note(sprintf("Dotplot %s: top %d term gösterildi.", out_base, nrow(df)))
}

plot_enrichment_dotplot(ORA_GO_FC_FILE,
                        "ora_go_fc_dotplot_publication",
                        "ORA — GO Biological Process (Frontal Cortex)",
                        mode = "ora")

plot_enrichment_dotplot(GSEA_GO_FC_FILE,
                        "gsea_go_fc_dotplot_publication",
                        "GSEA — GO Biological Process (Frontal Cortex)",
                        mode = "gsea")

plot_enrichment_dotplot(GSEA_KEGG_FC_FILE,
                        "gsea_kegg_fc_dotplot_publication",
                        "GSEA — KEGG Pathways (Frontal Cortex)",
                        mode = "gsea")


# -----------------------------------------------------------------------------
# 9. FIGURE NOTES
# -----------------------------------------------------------------------------
fig_notes <- c(
  "Figure notes — publication-ready outputs",
  paste(rep("=", 70), collapse = ""),
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("DESeq2 version: ", deseq2_version),
  paste0("Significance cutoffs: padj < ", PADJ_CUTOFF,
         " ; |log2FC| > ", LFC_CUTOFF),
  "",
  "Volcano plots (volcano_fc_publication.*, volcano_cb_publication.*)",
  paste0("  Input:  ", RES_FC_FILE, " / ", RES_CB_FILE),
  "  X-axis: log2FoldChange",
  "  Y-axis: -log10(adjusted p-value)",
  "  Horizontal dashed line: -log10(padj cutoff)",
  "  Vertical dashed lines:  +/- log2FC cutoff",
  "  Genes coloured by regulation status; top genes labelled with ggrepel.",
  "",
  "PCA plot (pca_plot_publication.*)",
  paste0("  Input:  ", DDS_FILE),
  "  Transformation: vst(dds, blind = FALSE)",
  "  Aesthetic: colour = diagnosis, shape = brain_region",
  "  Axes annotated with explained variance percentage.",
  "",
  "Heatmap (heatmap_top_genes_publication.*)",
  paste0("  Input:  ", DDS_FILE, " + ", RES_FC_FILE, " + ", METADATA_FILE),
  paste0("  Gene selection: top ", N_HEATMAP,
         " genes by ascending padj from the FC contrast"),
  "  Values:   vst counts, row-wise z-score scaled",
  "  Columns:  ordered (Control -> AD), then by brain_region; clustering off",
  "  Rows:     hierarchically clustered",
  "",
  "Enrichment dotplots",
  paste0("  Inputs: ", ORA_GO_FC_FILE, ", ", GSEA_GO_FC_FILE,
         ", ", GSEA_KEGG_FC_FILE),
  paste0("  Up to top ", N_DOTPLOT, " terms by ascending p.adjust"),
  "  Pathway names wrapped with stringr::str_wrap (width = 45)",
  "  Style consistent across ORA / GSEA panels.",
  "",
  "Improvements over previous figures",
  "  - Threshold lines drawn on volcano plots (padj and |log2FC|).",
  "  - Up / Down / Not significant colour coding standardised.",
  "  - Gene labels added via ggrepel (no overlap, top significant only).",
  "  - PCA axes display percent variance explained.",
  "  - Heatmap samples grouped by diagnosis (Control vs AD).",
  "  - Dotplot pathway names wrapped, top-N restricted, font size tuned.",
  "  - All outputs saved as 300 dpi PNG and vector PDF.",
  "",
  "Per-figure run-time notes:"
)
if (length(.notes) > 0) {
  fig_notes <- c(fig_notes, paste0("  - ", .notes))
}
writeLines(fig_notes, file.path(VIS_DIR, "figure_notes.txt"))
mark_done(file.path(VIS_DIR, "figure_notes.txt"))


# -----------------------------------------------------------------------------
# 10. ÖZET
# -----------------------------------------------------------------------------
cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
cat(" figure_generation_script.R — completed\n")
cat(paste(rep("=", 70), collapse = ""), "\n", sep = "")
cat(sprintf("DESeq2 version:        %s\n", deseq2_version))
cat(sprintf("Thresholds:            padj < %g  |  |log2FC| > %g\n",
            PADJ_CUTOFF, LFC_CUTOFF))
cat(sprintf("Labelled top genes:    %d (volcano)\n", N_LABEL_GENES))
cat(sprintf("Heatmap top genes:     %d\n",            N_HEATMAP))
cat(sprintf("Dotplot top terms:     %d\n",            N_DOTPLOT))
cat("\nProduced files:\n")
for (f in .produced) cat("  ", f, "\n", sep = "")
if (length(.skipped) > 0) {
  cat("\nSkipped figures:\n")
  for (s in .skipped) cat("  ", s, "\n", sep = "")
}
cat("\nNotes written to:\n")
cat("  ", file.path(VIS_DIR, "deseq2_version.txt"), "\n", sep = "")
cat("  ", file.path(VIS_DIR, "methods_notes.txt"),  "\n", sep = "")
cat("  ", file.path(VIS_DIR, "figure_notes.txt"),   "\n", sep = "")
cat(paste(rep("=", 70), collapse = ""), "\n\n", sep = "")