# for 
#     wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/977/SRR14513977/SRR14513977_1.fastq.gz
#     wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/977/SRR14513977/SRR14513977_2.fastq.gz
#     ./pipeline.sh fastq/SRR14513977
#     rm -rf fastq/SRR14513977_*

#!/bin/bash

TOPLAM=$(wc -l < wget_commands.txt)
SRR_SAYISI=$((TOPLAM / 2))

echo "=========================================="
echo "  Toplam: $SRR_SAYISI örnek işlenecek"
echo "  Her biri için:"
echo "    1) _1 ve _2 fastq indirilecek"
echo "    2) pipeline.sh çalıştırılacak"
echo "    3) fastq dosyaları silinecek"
echo "=========================================="
read -p "Başlamak istiyor musunuz? (e/h): " ONAY

if [[ "$ONAY" != "e" && "$ONAY" != "E" ]]; then
    echo "İptal edildi."
    exit 0
fi

echo ""
mkdir -p fastq

# wget_commands.txt içindeki satırları ikişer ikişer oku
while read line1 && read line2; do
    # _1.fastq.gz URL'sinden SRR adını çıkar
    SRR=$(basename $(echo "$line1" | awk '{print $NF}') _1.fastq.gz)

    echo "=========================================="
    echo "İndiriliyor: $SRR"
    echo "=========================================="

    # İndirme hedefini fastq/ klasörü olarak ayarla
    # wget -nc -P fastq/ $(echo "$line1" | awk '{print $NF}')
    # wget -nc -P fastq/ $(echo "$line2" | awk '{print $NF}')

    echo "Pipeline çalıştırılıyor: $SRR"
    ./pipeline.sh fastq/$SRR

    exit 0 # test için erken çıkış

    echo "Temizleniyor: fastq/${SRR}_*"
    rm -f fastq/${SRR}_1.fastq.gz fastq/${SRR}_2.fastq.gz


    echo "Tamamlandı: $SRR"
    echo ""

    # sleep inf # Test için sonsuz bekleme

done < wget_commands.txt