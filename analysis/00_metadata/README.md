cat > 00_metadata/README.md <<'EOF'
# 00_metadata

## Amaç
Bu adımda downstream RNA-seq analizi için gerekli örnek metadata tablosu hazırlanır.

## Girdi
- ../SraRunTable_clean.csv

## Kullanılan sütunlar
- Run
- AGE
- sex
- pmi
- brain_region
- diagnosis

## Yapılan işlemler
1. Analiz için gerekli sütunlar seçildi.
2. `Run` sütunu `sample` olarak yeniden adlandırıldı.
3. `AGE` sütunu `age` olarak yeniden adlandırıldı.
4. Sayısal sütunlar uygun formata çevrildi.
5. Metin alanlarında boşluk temizliği yapıldı.
6. Tekrarlı sample kayıtları kontrol edildi.
7. DESeq2 için temiz metadata dosyası üretildi.

## Çıktılar
- metadata_raw_selected.csv
- metadata.csv
- metadata_qc.txt

## Not
`metadata.csv` içindeki `sample` isimleri count matrix kolon isimleri ile birebir eşleşmelidir.
EOF



### 00_metadata için sonuç yorumu

CB: 96, FC: 95
Control: 96, AD: 95
male: 100, female: 91



Metadata QC Summary
===================

Total rows in original table: 191
Total unique samples: 191
Number of duplicated samples: 0

Missing values per column:
sample          0
age             0
sex             0
pmi             4
brain_region    0
diagnosis       0

Value counts for sex:
sex
male      100
female     91

Value counts for brain_region:
brain_region
CB    96
FC    95

Value counts for diagnosis:
diagnosis
Control    96
AD         95
