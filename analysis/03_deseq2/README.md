# Differential Expression Analysis (DESeq2)

## Amaç
Gen ekspresyon farklılıklarını analiz etmek.

## Model
~ AGE + sex + pmi + brain_region + diagnosis + brain_region:diagnosis

## Açıklama
- AGE, sex, pmi → confounder kontrolü
- brain_region → FC vs CB farkı
- diagnosis → AD vs Control farkı
- interaction → AD etkisi bölgeye göre değişiyor mu?

## Analizler
- AD vs Control (FC içinde)
- AD vs Control (CB içinde)
- Region difference
- Interaction effect

## Output
- DESeq2_results.csv
- significant_DEGs.csv