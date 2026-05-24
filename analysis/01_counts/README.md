# Count Matrix Preparation

## Amaç
featureCounts çıktılarından tek bir birleşik count matrix oluşturmak.

## Input
- ../counts/*_counts.txt

## Output
- merged_counts.tsv

## Açıklama
Her sample için ayrı oluşturulan count dosyaları birleştirilerek:
- satır = gen
- sütun = sample olacak şekilde tek tablo oluşturulur

## Not
Bu dosya DESeq2 için ana girdidir.