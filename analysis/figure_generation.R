#!/usr/bin/env Rscript
# =============================================================================
# figure_generation_script_poster.R
# -----------------------------------------------------------------------------
# Purpose
#   Bitirme projesi POSTERİ (90 cm x 70 cm) için optimize edilmiş, vektörel
#   (SVG) figürler üretir. Mevcut DESeq2 / enrichment CSV çıktılarını okur,
#   raw data ve tablolara DOKUNMAZ; sadece yeni klasöre yeni figürler yazar.
#
#   Eski 'figure_generation_script.R' dosyasından farkları:
#     - Tüm çıktılar 'yeni_figur/' klasörüne yazılır.
#     - Format: SVG (vektörel, posterde büyütünce piksellenmez)
#                + 300 dpi PNG yedeği.
#     - Yazı boyutları posterde 1-2 m mesafeden okunacak şekilde büyütüldü
#       (base_size = 24, eksen etiketleri ~22 pt, gen/term yazıları ~16-20 pt).
#     - Label çakışmaları için ggrepel parametreleri (force, box.padding,
#       max.overlaps) sıkılaştırıldı; gen isimleri kalın yazıldı.
#     - Dotplot pathway isimleri için str_wrap genişliği daraltıldı,
#       satır yüksekliği büyütüldü.
#     - Heatmap satır (gen) fontu 7 -> 14'e çıkarıldı; sütun isimleri 100+
#       örnek olduğunda zaten okunamayacağı için gizlendi, annotation barı
#       grup bilgisini taşıyor.
#     - Figür boyutları posterdeki yerleşime göre büyütüldü.
#
# Çalıştırma
#   cd ~/Pipeline_data_analysis/analysis
#   Rscript figure_generation_script_poster.R
# =============================================================================


# -----------------------------------------------------------------------------
# 0. KONFİGÜRASYON
# -----------------------------------------------------------------------------
# Yollar (girdiler eskisi gibi)
METADATA_FILE      <- "00_metadata/metadata_filtered.csv"
DDS_FILE           <- "03_deseq2/dds.rds"
RES_FC_FILE        <- "03_deseq2/results_diagnosis_in_FC.csv"
RES_CB_FILE        <- "03_deseq2/results_diagnosis_in_CB.csv"

ORA_GO_FC_FILE     <- "05_enrichment/ora_go_fc.csv"
GSEA_GO_FC_FILE    <- "05_enrichment/gsea_go_fc.csv"
GSEA_KEGG_FC_FILE  <- "05_enrichment/gsea_kegg_fc.csv"

# >>> YENİ ÇIKTI KLASÖRÜ <<<
OUT_DIR <- "yeni_figur"

# Significance thresholds (tek noktada)
PADJ_CUTOFF   <- 0.05
LFC_CUTOFF    <- 1
N_LABEL_GENES <- 12          # volcano'da etiketlenecek gen sayısı
N_HEATMAP     <- 50          # heatmap'te gösterilecek top gen sayısı
N_DOTPLOT     <- 15          # dotplot'ta gösterilecek top pathway sayısı

# Poster için temel font boyutu (tüm tema buradan ölçeklenir)
POSTER_BASE_SIZE <- 24

# Üretilen / atlanan görsel kayıtları
.produced <- character(0)
.skipped  <- character(0)
.notes    <- character(0)
add_note  <- function(msg) .notes <<- c(.notes, msg)
mark_done <- function(...) .produced <<- c(.produced, c(...))
mark_skip <- function(name, reason) .skipped <<- c(.skipped, paste0(name, "  --  ", reason))

# Çıktı klasörü
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)


# -----------------------------------------------------------------------------
# 1. PAKET KONTROLÜ
# -----------------------------------------------------------------------------
required_pkgs <- c("DESeq2", "ggplot2", "ggrepel", "pheatmap",
                   "dplyr", "readr", "stringr", "tibble")

missing_pkgs <- required_pkgs[!vapply(required_pkgs,
                                      requireNamespace, logical(1),
                                      quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Eksik paket(ler): ", paste(missing_pkgs, collapse = ", "), "\n",
    "Kurulum:\n",
    "  CRAN:        install.packages(c('ggplot2','ggrepel','pheatmap',",
    "'dplyr','readr','stringr','tibble','svglite'))\n",
    "  Bioconductor: BiocManager::install(c('DESeq2'))\n"
  )
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

# Opsiyonel: ENSEMBL -> SYMBOL haritası
has_orgdb <- requireNamespace("org.Hs.eg.db", quietly = TRUE) &&
             requireNamespace("AnnotationDbi", quietly = TRUE)

# Opsiyonel: svglite (varsa font handling'i daha iyi yapar, özellikle TR karakter)
has_svglite <- requireNamespace("svglite", quietly = TRUE)


# -----------------------------------------------------------------------------
# 2. POSTER TEMASI (büyük fontlar, kalın başlıklar, ferah margin)
# -----------------------------------------------------------------------------
theme_poster <- function(base_size = POSTER_BASE_SIZE) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major   = element_line(colour = "grey92", linewidth = 0.4),
      plot.title = element_text(face = "bold",
                                hjust = 0.5,
                                size = base_size + 2,
                                margin = margin(b = 14)),
      plot.subtitle      = element_text(hjust = 0.5,
                                        size = base_size + 1,
                                        margin = margin(b = 12)),
      axis.title         = element_text(face = "bold", size = base_size + 2),
      axis.title.x       = element_text(margin = margin(t = 12)),
      axis.title.y       = element_text(margin = margin(r = 12)),
      axis.text          = element_text(colour = "black", size = base_size),
      axis.ticks         = element_line(linewidth = 0.5),
      legend.title       = element_text(face = "bold", size = base_size + 1),
      legend.text        = element_text(size = base_size),
      legend.position    = "right",
      legend.key.size    = unit(1.4, "lines"),
      legend.background  = element_blank(),
      strip.background   = element_rect(fill = "grey95", colour = NA),
      strip.text         = element_text(face = "bold", size = base_size + 1),
      plot.margin = margin(28, 35, 28, 35)
    )
}

# Tutarlı renkler
COL_UP    <- "#C0392B"   # kırmızı  - upregulated
COL_DOWN  <- "#1F618D"   # mavi     - downregulated
COL_NS    <- "grey70"    # gri      - not significant
COL_GROUP <- c("Control" = "#2E86AB", "AD" = "#E63946")


# -----------------------------------------------------------------------------
# 3. KAYIT FONKSİYONU — SVG (vektörel) + PNG (yedek 300 dpi)
# -----------------------------------------------------------------------------
save_pub <- function(plot_obj, file_base, width = 12, height = 9) {
  svg_path <- file.path(OUT_DIR, paste0(file_base, ".svg"))
  png_path <- file.path(OUT_DIR, paste0(file_base, ".png"))

  # SVG (vektörel — posterde büyütünce kalite kaybı YOK)
  if (has_svglite) {
    ggsave(svg_path, plot_obj, device = svglite::svglite,
           width = width, height = height, units = "in")
  } else {
    ggsave(svg_path, plot_obj, device = "svg",
           width = width, height = height, units = "in")
  }

  # PNG yedek (300 dpi — print için bir tık ekstra güven)
  ggsave(png_path, plot_obj,
         width = width, height = height, dpi = 300, units = "in")

  mark_done(svg_path, png_path)
}


# -----------------------------------------------------------------------------
# 4. YARDIMCI FONKSİYONLAR
# -----------------------------------------------------------------------------
prepare_de_table <- function(df, source_label) {
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
  df$label <- ifelse(is.na(df$symbol) | df$symbol == "",
                     df$ensembl_clean, df$symbol)
  df
}


# -----------------------------------------------------------------------------
# 5. VOLCANO PLOT (poster-ready)
# -----------------------------------------------------------------------------
make_volcano <- function(res_path, out_base, title_text) {
  if (!file.exists(res_path)) {
    mark_skip(out_base, sprintf("Girdi dosyası bulunamadı: %s", res_path))
    return(invisible(NULL))
  }

  df <- read.csv(res_path, stringsAsFactors = FALSE, check.names = FALSE)
  df <- prepare_de_table(df, basename(res_path))
  df <- add_gene_symbols(df)

  df <- df[!is.na(df$padj) & !is.na(df$log2FoldChange), ]

  df$category <- "Not significant"
  df$category[df$padj < PADJ_CUTOFF & df$log2FoldChange >  LFC_CUTOFF] <- "Upregulated"
  df$category[df$padj < PADJ_CUTOFF & df$log2FoldChange < -LFC_CUTOFF] <- "Downregulated"
  df$category <- factor(df$category,
                        levels = c("Upregulated", "Downregulated", "Not significant"))

  sig_df <- df[df$category != "Not significant", ]
  sig_df <- sig_df[order(sig_df$padj), ]
  top_df <- head(sig_df, N_LABEL_GENES)

  max_neglog <- max(-log10(df$padj[df$padj > 0]), na.rm = TRUE)
  cap_val <- max_neglog * 1.05
  df$neglog10padj     <- pmin(-log10(df$padj), cap_val)
  top_df$neglog10padj <- pmin(-log10(top_df$padj), cap_val)

  p <- ggplot(df, aes(x = log2FoldChange, y = neglog10padj, colour = category)) +
    geom_point(alpha = 0.65, size = 2.6) +                       # büyük noktalar
    geom_hline(yintercept = -log10(PADJ_CUTOFF),
               linetype = "dashed", colour = "grey30", linewidth = 0.6) +
    geom_vline(xintercept = c(-LFC_CUTOFF, LFC_CUTOFF),
               linetype = "dashed", colour = "grey30", linewidth = 0.6) +
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
    theme_poster() +
    # Legend noktalarını da büyüt (yoksa küçücük kalıyor)
    guides(colour = guide_legend(override.aes = list(size = 6, alpha = 1)))

  if (nrow(top_df) > 0) {
    p <- p +
      geom_text_repel(
        data            = top_df,
        aes(label = label),
        size            = 7,                # ~18pt — uzaktan okunur
        fontface        = "bold",
        colour          = "black",
        bg.colour       = "white",          # yazıyı beyaz halka ile çevreler
        bg.r            = 0.15,
        box.padding     = 0.9,
        point.padding   = 0.5,
        segment.colour  = "grey25",
        segment.size    = 0.5,
        max.overlaps    = Inf,
        min.segment.length = 0,
        force           = 4,
        force_pull      = 0.4,
        seed            = 42,
        show.legend     = FALSE
      )
  }

  save_pub(p, out_base, width = 12, height = 9.5)

  add_note(sprintf(
    "Volcano (%s): %d up, %d down, %d not significant; %d genes labelled.",
    out_base,
    sum(df$category == "Upregulated"),
    sum(df$category == "Downregulated"),
    sum(df$category == "Not significant"),
    nrow(top_df)
  ))
}

make_volcano(RES_FC_FILE, "volcano_fc_poster",
             "Volcano plot — AD vs. Control in Frontal Cortex (FC)")
make_volcano(RES_CB_FILE, "volcano_cb_poster",
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
    mark_skip("pca_plot_poster",
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
      geom_point(size = 5, alpha = 0.85, stroke = 0.8) +        # büyük noktalar
      scale_colour_manual(values = COL_GROUP, name = "Diagnosis") +
      scale_shape_manual(values = c("CB" = 16, "FC" = 17),
                         name   = "Brain region") +
      labs(
        title = "PCA Plot of Variance-Stabilized Expression Data",
        x     = sprintf("PC1 (%.1f%% variance)", pct_var[1]),
        y     = sprintf("PC2 (%.1f%% variance)", pct_var[2])
      ) +
      theme_poster() +
      guides(colour = guide_legend(order = 1,
                                   override.aes = list(size = 6)),
             shape  = guide_legend(order = 2,
                                   override.aes = list(size = 6)))

    save_pub(p_pca, "pca_plot_poster", width = 12, height = 9)
    add_note(sprintf(
      "PCA: vst(dds, blind = FALSE); %d sample, PC1 = %.1f%%, PC2 = %.1f%%.",
      nrow(pca_data), pct_var[1], pct_var[2]
    ))
  }
} else {
  mark_skip("pca_plot_poster", "dds.rds yok.")
}


# -----------------------------------------------------------------------------
# 7. HEATMAP (poster-ready)
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

  res_fc <- res_fc[!is.na(res_fc$padj), ]
  res_fc <- res_fc[order(res_fc$padj), ]
  top_genes <- res_fc$gene_id[seq_len(min(N_HEATMAP, nrow(res_fc)))]

  vst_mat   <- assay(vsd)
  top_genes <- intersect(top_genes, rownames(vst_mat))

  if (length(top_genes) >= 2) {
    mat <- vst_mat[top_genes, , drop = FALSE]
    mat <- t(scale(t(mat)))                  # row-wise z-score

    # Satır isimlerini SYMBOL'a çevir
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
        new_names <- make.unique(new_names)
        rownames(mat) <- new_names
      }
    }

    ann_cols <- intersect(c("diagnosis", "brain_region"), colnames(meta_full))
    ann_df   <- meta_full[colnames(mat), ann_cols, drop = FALSE]
    ann_df$diagnosis <- factor(ann_df$diagnosis,
                               levels = c("Control", "AD"))
    if ("brain_region" %in% colnames(ann_df)) {
      ann_df$brain_region <- factor(ann_df$brain_region,
                                    levels = c("CB", "FC"))
    }

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

    # Posterde 100+ örnek ID'sini okumak imkansız - sütun isimlerini gizle,
    # annotation barı (diagnosis + brain_region) grupları zaten gösteriyor.
    show_colnames <- (ncol(mat) <= 20)

    draw_hm <- function() {
      pheatmap(
        mat,
        annotation_col    = ann_df,
        annotation_colors = ann_colors,
        cluster_cols      = FALSE,
        cluster_rows      = TRUE,
        show_colnames     = show_colnames,
        show_rownames     = TRUE,
        fontsize          = 16,           # genel font (legend, başlık)
        fontsize_row      = 14,           # gen isimleri — okunur büyüklük
        fontsize_col      = 11,
        treeheight_row    = 35,
        treeheight_col    = 0,
        border_color      = NA,
        color             = colorRampPalette(
                              c("#1F618D", "white", "#C0392B"))(100),
        annotation_legend = TRUE,
        legend            = TRUE,
        main              = sprintf(
          "Top %d DE genes (FC contrast) — vst, row z-score",
          length(top_genes)),
        silent            = TRUE
      )
    }

    hm_obj <- draw_hm()

    hm_svg <- file.path(OUT_DIR, "heatmap_top_genes_poster.svg")
    hm_png <- file.path(OUT_DIR, "heatmap_top_genes_poster.png")

    # SVG (vektörel)
    if (has_svglite) {
      svglite::svglite(hm_svg, width = 14, height = 17)
    } else {
      svg(hm_svg, width = 14, height = 17)
    }
    grid::grid.newpage(); grid::grid.draw(hm_obj$gtable)
    dev.off()

    # PNG yedek
    png(hm_png, width = 14, height = 17, units = "in", res = 300)
    grid::grid.newpage(); grid::grid.draw(hm_obj$gtable)
    dev.off()

    mark_done(hm_svg, hm_png)
    add_note(sprintf(
      "Heatmap: top %d genes by FC padj; column clustering off; samples ordered Control -> AD (by brain_region).",
      length(top_genes)
    ))
  } else {
    mark_skip("heatmap_top_genes_poster",
              "vsd'de eşleşen top gene sayısı < 2.")
  }
} else {
  mark_skip("heatmap_top_genes_poster",
            "Gerekli dosyalar eksik (dds.rds / FC results / metadata).")
}


# -----------------------------------------------------------------------------
# 8. ENRICHMENT DOTPLOTLARI (ORA GO, GSEA GO, GSEA KEGG) — poster-ready
# -----------------------------------------------------------------------------
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
    x_label    <- "Gene ratio"
    size_label <- "Count"
  } else {
    df$x_axis_value <- if ("NES" %in% colnames(df)) df$NES else NA_real_
    df$size_value   <- if ("setSize" %in% colnames(df)) df$setSize else NA_real_
    x_label    <- "Normalized enrichment score (NES)"
    size_label <- "Set size"
  }

  df <- df[!is.na(df$p.adjust), ]
  df <- df[order(df$p.adjust), ]
  df <- head(df, N_DOTPLOT)

  if (nrow(df) == 0) {
    mark_skip(out_base, "Çizilecek pathway kalmadı.")
    return(invisible(NULL))
  }

  # Pathway isimleri 40 karakterden uzunsa sar (poster fontunda hala okunur)
  df$Description_wrap <- str_wrap(df$Description, width = 40)
  df$Description_wrap <- factor(df$Description_wrap,
                                levels = rev(df$Description_wrap))

  p <- ggplot(df, aes(x = x_axis_value, y = Description_wrap)) +
    geom_point(aes(size = size_value, colour = p.adjust)) +
    scale_colour_gradient(low = "#C0392B", high = "#1F618D",
                          name = "p.adjust",
                          guide = guide_colourbar(reverse = FALSE,
                                                  barwidth  = 1.4,
                                                  barheight = 14)) +
    scale_size_continuous(name = size_label, range = c(5, 16)) +
    labs(title = title_text, x = x_label, y = NULL) +
    theme_poster(base_size = 22) +
    theme(
      # Pathway adları iki satıra sığabiliyor — lineheight ve size öncelikli
      axis.text.y = element_text(size = 18, lineheight = 0.95,
                                 colour = "black"),
      panel.grid.major.y = element_line(colour = "grey94", linewidth = 0.3)
    )

  if (mode == "gsea") {
    p <- p + geom_vline(xintercept = 0, linetype = "dashed",
                        colour = "grey40", linewidth = 0.5)
  }

  # Boy: 15 pathway için ~15 inch civarı (her terim ~0.9-1 inch)
  fig_h <- 6.5 + 0.55 * nrow(df)
  save_pub(p, out_base, width = 14, height = fig_h)

  add_note(sprintf("Dotplot %s: top %d term gösterildi.", out_base, nrow(df)))
}

plot_enrichment_dotplot(ORA_GO_FC_FILE,
                        "ora_go_fc_dotplot_poster",
                        "ORA — GO Biological Process (Frontal Cortex)",
                        mode = "ora")

plot_enrichment_dotplot(GSEA_GO_FC_FILE,
                        "gsea_go_fc_dotplot_poster",
                        "GSEA — GO Biological Process (Frontal Cortex)",
                        mode = "gsea")

plot_enrichment_dotplot(GSEA_KEGG_FC_FILE,
                        "gsea_kegg_fc_dotplot_poster",
                        "GSEA — KEGG Pathways (Frontal Cortex)",
                        mode = "gsea")


# -----------------------------------------------------------------------------
# 9. NOT DOSYALARI
# -----------------------------------------------------------------------------
deseq2_version <- as.character(packageVersion("DESeq2"))

fig_notes <- c(
  "Figure notes — POSTER versiyonu (90 cm x 70 cm)",
  paste(rep("=", 70), collapse = ""),
  "",
  paste0("Generated:      ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("DESeq2 version: ", deseq2_version),
  paste0("R version:      ", paste(R.version$major, R.version$minor, sep = ".")),
  paste0("svglite used:   ", has_svglite),
  paste0("org.Hs.eg.db:   ", has_orgdb),
  "",
  "Output folder:  ", paste0("  ", OUT_DIR),
  "",
  "Format",
  "  SVG (vektörel — posterde büyütünce piksellenmez)",
  "  PNG 300 dpi (yedek; print için güvence)",
  "",
  "Tipografi",
  paste0("  Base font size: ", POSTER_BASE_SIZE, " pt"),
  "  Volcano label fontu: ~18 pt, bold, beyaz halo",
  "  Heatmap satır (gen) fontu: 14 pt",
  "  Dotplot pathway fontu:     18 pt, str_wrap(width = 40)",
  "",
  "Threshold'lar",
  paste0("  padj < ", PADJ_CUTOFF, " ; |log2FC| > ", LFC_CUTOFF),
  paste0("  Volcano label sayısı: ", N_LABEL_GENES),
  paste0("  Heatmap top gene:     ", N_HEATMAP),
  paste0("  Dotplot top term:     ", N_DOTPLOT),
  "",
  "Per-figür kayıtlar:"
)
if (length(.notes) > 0) {
  fig_notes <- c(fig_notes, paste0("  - ", .notes))
}
writeLines(fig_notes, file.path(OUT_DIR, "figure_notes_poster.txt"))
mark_done(file.path(OUT_DIR, "figure_notes_poster.txt"))


# -----------------------------------------------------------------------------
# 10. ÖZET
# -----------------------------------------------------------------------------
cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
cat(" figure_generation_script_poster.R — tamamlandı\n")
cat(paste(rep("=", 70), collapse = ""), "\n", sep = "")
cat(sprintf("Output klasörü:        %s/\n", OUT_DIR))
cat(sprintf("Format:                SVG (vektörel) + PNG 300 dpi yedek\n"))
cat(sprintf("DESeq2 version:        %s\n", deseq2_version))
cat(sprintf("svglite kullanıldı:    %s\n", has_svglite))
cat(sprintf("Thresholds:            padj < %g  |  |log2FC| > %g\n",
            PADJ_CUTOFF, LFC_CUTOFF))
cat("\nÜretilen dosyalar:\n")
for (f in .produced) cat("  ", f, "\n", sep = "")
if (length(.skipped) > 0) {
  cat("\nAtlanan figürler:\n")
  for (s in .skipped) cat("  ", s, "\n", sep = "")
}
cat(paste(rep("=", 70), collapse = ""), "\n\n", sep = "")