# # # for 
# # #     wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/977/SRR14513977/SRR14513977_1.fastq.gz
# # #     wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/977/SRR14513977/SRR14513977_2.fastq.gz
# # #     ./pipeline.sh fastq/SRR14513977
# # #     rm -rf fastq/SRR14513977_*

# # #!/bin/bash

# # TOPLAM=$(wc -l < wget_commands.txt)
# # SRR_SAYISI=$((TOPLAM / 2))

# # echo "=========================================="
# # echo "  Toplam: $SRR_SAYISI örnek işlenecek"
# # echo "  Her biri için:"
# # echo "    1) _1 ve _2 fastq indirilecek"
# # echo "    2) pipeline.sh çalıştırılacak"
# # echo "    3) fastq dosyaları silinecek"
# # echo "=========================================="
# # read -p "Başlamak istiyor musunuz? (e/h): " ONAY

# # if [[ "$ONAY" != "e" && "$ONAY" != "E" ]]; then
# #     echo "İptal edildi."
# #     exit 0
# # fi

# # echo ""
# # mkdir -p fastq

# # # wget_commands.txt içindeki satırları ikişer ikişer oku
# # while read line1 && read line2; do
# #     # _1.fastq.gz URL'sinden SRR adını çıkar
# #     SRR=$(basename $(echo "$line1" | awk '{print $NF}') _1.fastq.gz)

# #     echo "=========================================="
# #     echo "İndiriliyor: $SRR"
# #     echo "=========================================="

# #     # İndirme hedefini fastq/ klasörü olarak ayarla
# #     wget -nc -P fastq/ $(echo "$line1" | awk '{print $NF}')
# #     wget -nc -P fastq/ $(echo "$line2" | awk '{print $NF}')

# #     echo "Pipeline çalıştırılıyor: $SRR"
# #     ./pipeline.sh $SRR

# #     exit 0 # test için erken çıkış

# #     echo "Temizleniyor: fastq/${SRR}_*"
# #     rm -f fastq/${SRR}_1.fastq.gz fastq/${SRR}_2.fastq.gz


# #     echo "Tamamlandı: $SRR"
# #     echo ""

# #     # sleep inf # Test için sonsuz bekleme

# # done < wget_commands.txt


# #!/bin/bash
# set -e
# set -o pipefail

# TOPLAM=$(grep -c '^wget ' wget_commands.txt)
# SRR_SAYISI=$((TOPLAM / 2))

# echo "=========================================="
# echo "  Toplam: $SRR_SAYISI örnek işlenecek"
# echo "=========================================="
# read -p "Başlamak istiyor musunuz? (e/h): " ONAY

# if [[ "$ONAY" != "e" && "$ONAY" != "E" ]]; then
#     echo "İptal edildi."
#     exit 0
# fi

# mkdir -p fastq

# while read -r line1 && read -r line2; do
#     url1=$(echo "$line1" | awk '{print $NF}')
#     url2=$(echo "$line2" | awk '{print $NF}')
#     SRR=$(basename "$url1" _1.fastq.gz)

#     echo "=========================================="
#     echo "İndiriliyor: $SRR"
#     echo "=========================================="

#     wget -nv -nc -P fastq/ "$url1"
#     wget -nv -nc -P fastq/ "$url2"

#     ls -lh "fastq/${SRR}_1.fastq.gz" "fastq/${SRR}_2.fastq.gz"

#     echo "Pipeline çalıştırılıyor: $SRR"
#     ./pipeline.sh "$SRR"

#     # exit 0   # test için bilerek kalsın

# done < <(grep '^wget ' wget_commands.txt)


#!/bin/bash
set -euo pipefail

WGET_FILE="wget_commands.txt"
PIPELINE_SCRIPT="./pipeline.sh"
FASTQ_DIR="fastq"

if [ ! -f "$WGET_FILE" ]; then
  echo "ERROR: $WGET_FILE not found"
  exit 1
fi

if [ ! -f "$PIPELINE_SCRIPT" ]; then
  echo "ERROR: $PIPELINE_SCRIPT not found"
  exit 1
fi

if [ ! -x "$PIPELINE_SCRIPT" ]; then
  chmod +x "$PIPELINE_SCRIPT"
fi

TOPLAM=$(grep -c '^wget ' "$WGET_FILE" || true)

if [ "$TOPLAM" -eq 0 ]; then
  echo "ERROR: No wget lines found in $WGET_FILE"
  exit 1
fi

if [ $((TOPLAM % 2)) -ne 0 ]; then
  echo "ERROR: wget_commands.txt contains an odd number of wget lines."
  echo "Each sample must have exactly 2 lines (_1 and _2)."
  exit 1
fi

SRR_SAYISI=$((TOPLAM / 2))

echo "=========================================="
echo "  Toplam: $SRR_SAYISI örnek işlenecek"
echo "  Wget dosyası: $WGET_FILE"
echo "=========================================="
read -r -p "Başlamak istiyor musunuz? (e/h): " ONAY

if [[ "$ONAY" != "e" && "$ONAY" != "E" ]]; then
  echo "İptal edildi."
  exit 0
fi

mkdir -p "$FASTQ_DIR"

while read -r line1 && read -r line2; do
  url1=$(echo "$line1" | awk '{print $NF}')
  url2=$(echo "$line2" | awk '{print $NF}')

  SRR=$(basename "$url1" _1.fastq.gz)

  if [[ -z "$SRR" ]]; then
    echo "ERROR: Could not parse SRR from line: $line1"
    exit 1
  fi

  echo "=========================================="
  echo "İndiriliyor: $SRR"
  echo "=========================================="

  wget -nv -nc -P "$FASTQ_DIR/" "$url1"
  wget -nv -nc -P "$FASTQ_DIR/" "$url2"

  if [ ! -f "$FASTQ_DIR/${SRR}_1.fastq.gz" ] || [ ! -f "$FASTQ_DIR/${SRR}_2.fastq.gz" ]; then
    echo "ERROR: FASTQ files missing after download for $SRR"
    exit 1
  fi

  ls -lh "$FASTQ_DIR/${SRR}_1.fastq.gz" "$FASTQ_DIR/${SRR}_2.fastq.gz"

  echo "Pipeline çalıştırılıyor: $SRR"
  "$PIPELINE_SCRIPT" "$SRR"

  echo "Temizleniyor: $FASTQ_DIR/${SRR}_*"
  rm -f \
    "$FASTQ_DIR/${SRR}_1.fastq.gz" \
    "$FASTQ_DIR/${SRR}_2.fastq.gz"

  echo "Tamamlandı: $SRR"
  echo

done < <(grep '^wget ' "$WGET_FILE")

echo "=========================================="
echo "Tüm sample'lar tamamlandı"
echo "=========================================="