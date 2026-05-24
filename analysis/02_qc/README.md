# 02_qc

## Amaç
Bu adımda sample kalite kontrolü yapılır ve düşük kaliteli örnekler downstream analizden çıkarılmak üzere belirlenir.

## Girdi
- ../counts/*_counts.txt.summary

## Yapılan işlemler
1. Tüm featureCounts summary dosyaları toplandı.
2. Her sample için assigned read sayısı hesaplandı.
3. Toplam read sayısı üzerinden assigned percent hesaplandı.
4. %35 threshold altında kalan sample'lar işaretlendi.
5. Analize dahil edilecek ve çıkarılacak örnek listeleri oluşturuldu.

## Çıktılar
- featurecounts_qc.tsv
- kept_samples.tsv
- removed_samples.tsv
- qc_summary.txt

## Not
Bu aşamada verilen %35 threshold, featureCounts assigned percent değerine göre uygulanır.

 ### Sonuç yorumu

Elindeki çıktıya göre:

toplam sample: 191
tutulan sample: 148
çıkarılan sample: 43

    66 AD
        22 CB
        44 FC
    82 Control
        36 CB
        46 FC


eksik pmıler çıkarıldı 

    65 AD
        22 CB
        43 FC
    80 Control
        35 CB
        45 FC