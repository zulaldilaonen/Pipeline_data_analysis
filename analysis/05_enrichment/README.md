# 05_enrichment

## Amaç
Bu adımda diferansiyel ekspresyon analizi sonrası fonksiyonel zenginleştirme analizi yapılır.

## Girdi
- ../03_deseq2/results_diagnosis_in_FC.csv

## Yapılan işlemler
1. FC için diferansiyel ekspresyon sonuçları okundu.
2. Ensembl gen kimlikleri temizlendi.
3. Anlamlı gen listesi çıkarıldı.
4. ORA analizi uygulandı:
   - GO Biological Process
   - KEGG
5. Tüm sıralı gen listesi ile GSEA uygulandı:
   - GO Biological Process
   - KEGG
6. Sonuç tabloları ve görseller kaydedildi.

## Çıktılar
- ora_go_fc.csv
- ora_kegg_fc.csv
- gsea_go_fc.csv
- gsea_kegg_fc.csv
- ora_go_fc_dotplot.png
- ora_kegg_fc_dotplot.png
- gsea_go_fc_dotplot.png
- gsea_kegg_fc_dotplot.png

## Not
Bu adım frontal cortex (FC) için yapılmaktadır. İstenirse daha sonra cerebellum (CB) için de aynı analiz uygulanabilir.



 # 1. ORA (Over-Representation Analysis)
Mantık:

“Anlamlı değişen genler hangi biyolojik süreçlerde beklenenden fazla yer alıyor?”

Senin scriptte:
```
padj < 0.05 ve |log2FC| > 1
```


→ sadece güçlü değişen genler kullanılıyor

 ## ora_go_fc.csv
- GO Biological Process (örn: inflammation, synapse, neuron death)
- Alzheimer’da hangi süreçler bozulmuş gösterir

 ## ora_kegg_fc.csv
- Pathway (örn: Alzheimer disease, oxidative phosphorylation)
- Hastalığın moleküler mekanizması


# 2. GSEA (Gene Set Enrichment Analysis)
Mantık:

“Tüm genleri sıralarsak, belirli biyolojik yolaklar yukarı mı aşağı mı kaymış?”

Farkı:

ORA → sadece seçilmiş genler
GSEA → tüm veri

## gsea_go_fc.csv
Süreçlerin genel aktivasyonu
## gsea_kegg_fc.csv
Pathway seviyesinde değişim


# 3. Önemli kolonlar (çok kritik)

CSV’lerde göreceğin:

Kolon	             Anlam

pvalue	         ham anlamlılık
p.adjust	         düzeltilmiş p (en önemlisi)
GeneRatio	      o pathway’deki gen oranı
NES (GSEA)	      yön + şiddet
Description	      süreç adı


## yorumlama 

### ORA için:
p.adjust < 0.05 → önemli
Örnek yorum:
"inflammatory response ↑"
"synaptic signaling ↓"
### GSEA için:
NES > 0 → upregulated
NES < 0 → downregulated