# Paketleri hangi dizine yüklüyor, ve hangi dizindeki paketler kullanılıyor
Rscript -e ".libPaths()"


# Paket analizi
echo "=== 03_deseq2/deseq2_analysis.R ===" && \
cat /home/talhaaygen/Pipeline_data_analysis/analysis/03_deseq2/deseq2_analysis.R && \
echo "=== 04_visualization/visualization.R ===" && \
cat /home/talhaaygen/Pipeline_data_analysis/analysis/04_visualization/visualization.R && \
echo "=== 05_enrichment/enrichment_fc.R ===" && \
cat /home/talhaaygen/Pipeline_data_analysis/analysis/05_enrichment/enrichment_fc.R

3 scripti analiz ettim. Tam paket listesi şu:

| Script | Paket | Kaynak |
|--------|-------|--------|
| deseq2_analysis.R | `DESeq2` | Bioconductor |
| deseq2_analysis.R | `apeglm` | Bioconductor |
| visualization.R | `DESeq2` | Bioconductor |
| visualization.R | `ggplot2` | CRAN |
| visualization.R | `pheatmap` | CRAN |
| enrichment_fc.R | `clusterProfiler` | Bioconductor |
| enrichment_fc.R | `org.Hs.eg.db` | Bioconductor |
| enrichment_fc.R | `enrichplot` | Bioconductor |
| enrichment_fc.R | `ggplot2` | CRAN |

`treeio` bug'ını da göz önünde bulundurarak Dockerfile'ı yazıyorum:Dockerfile hazır. Şimdi sunucuya kopyala ve build et:


# Hangi paketler kurulu? hangileri eksik?
Rscript -e "
pkgs <- c('DESeq2', 'apeglm', 'ggplot2', 'pheatmap', 'remotes',
          'treeio', 'ggtree', 'clusterProfiler', 'org.Hs.eg.db', 'enrichplot')
for (p in pkgs) {
  status <- requireNamespace(p, quietly=TRUE)
  cat(sprintf('%-20s %s\n', p, ifelse(status, 'OK', 'EKSIK')))
}
"


# Kurulum adımları

Rscript -e "install.packages(c('remotes', 'ggplot2', 'pheatmap'), repos='https://cloud.r-project.org', lib='/usr/local/lib/R/site-library')"

apt-get install -y libgit2-dev # bu paket zaten kuruluydu
Rscript -e "remotes::install_github('YuLab-SMU/treeio', upgrade='never', lib='/usr/local/lib/R/site-library')"

Rscript -e "remotes::install_github('YuLab-SMU/ggtree', upgrade='never', lib='/usr/local/lib/R/site-library')"

Rscript -e "BiocManager::install(c('DESeq2', 'apeglm'), ask=FALSE, update=FALSE)"

# gerekirse
Rscript -e 'install.packages("BH")'

# gerekirse
# 1. Flag'i ekle
echo 'CXXFLAGS = -O2 -fpermissive' >> ~/.R/Makevars

# 2. Kurulumu çalıştır
Rscript -e "BiocManager::install(c('clusterProfiler', 'org.Hs.eg.db', 'enrichplot'), ask=FALSE, update=FALSE, force=TRUE)"

# 3. Kurulum bitince flag'i temizle
sed -i '/fpermissive/d' ~/.R/Makevars


# ya da 
Rscript -e "remotes::install_version('BH', version='1.81.0-1')"
Rscript -e "install.packages('remotes')"

Rscript -e "BiocManager::install(c('clusterProfiler', 'org.Hs.eg.db', 'enrichplot'), ask=FALSE, update=FALSE)"


