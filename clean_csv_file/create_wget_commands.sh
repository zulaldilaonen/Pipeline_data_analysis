#!/bin/bash
set -euo pipefail

INPUT="SRR_Acc_List_clean.txt"
OUTPUT="wget_commands_clean.txt"

> "$OUTPUT"

while read -r srr; do
  [ -z "$srr" ] && continue

  prefix="${srr:0:6}"     # SRR145
  suffix="${srr: -3}"     # son 3 rakam, örn 977

  echo "wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${prefix}/${suffix}/${srr}/${srr}_1.fastq.gz" >> "$OUTPUT"
  echo "wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${prefix}/${suffix}/${srr}/${srr}_2.fastq.gz" >> "$OUTPUT"
done < "$INPUT"

echo "Oluşturuldu: $OUTPUT"
wc -l "$OUTPUT"
