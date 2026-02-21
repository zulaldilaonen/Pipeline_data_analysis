#!/bin/bash
set -e
set -o pipefail

########################################
# RNA-Seq Pipeline
# Usage: ./pipeline.sh SRRXXXXXXX
########################################

if [ -z "$1" ]; then
  echo "Usage: ./pipeline.sh SRRXXXXXXX"
  exit 1
fi

SAMPLE=$1
THREADS=8

BASE_DIR=$(pwd)

FASTQ_DIR="$BASE_DIR/fastq"
TRIM_DIR="$BASE_DIR/trimmed"
ALIGN_DIR="$BASE_DIR/alignment"
COUNT_DIR="$BASE_DIR/counts"
QC_DIR="$BASE_DIR/QC"
LOG_DIR="$BASE_DIR/logs"
REF_DIR="$BASE_DIR/reference"

GENOME_INDEX="$REF_DIR/GRCh38_index"
GTF="$REF_DIR/gencode.v45.annotation.gtf"
ADAPTERS="$REF_DIR/adapters.fa"
TRIMMOMATIC_JAR="/usr/share/java/trimmomatic.jar"

mkdir -p $FASTQ_DIR $TRIM_DIR $ALIGN_DIR $COUNT_DIR $QC_DIR $LOG_DIR $REF_DIR

echo "======================================"
echo "Processing sample: $SAMPLE"
echo "======================================"

########################################
# TOOL CHECK
########################################

for tool in fastqc multiqc hisat2 samtools featureCounts java; do
  if ! command -v $tool &> /dev/null; then
    echo "ERROR: $tool is not installed."
    exit 1
  fi
done

if [ ! -f "$TRIMMOMATIC_JAR" ]; then
  echo "ERROR: Trimmomatic jar not found at $TRIMMOMATIC_JAR"
  exit 1
fi

########################################
# INPUT CHECK
########################################

if [ ! -f "$FASTQ_DIR/${SAMPLE}_1.fastq.gz" ] || [ ! -f "$FASTQ_DIR/${SAMPLE}_2.fastq.gz" ]; then
  echo "ERROR: FASTQ files not found in $FASTQ_DIR"
  exit 1
fi

if [ ! -f "$GENOME_INDEX.1.ht2" ]; then
  echo "ERROR: HISAT2 index not found in reference/"
  exit 1
fi

########################################
# 1. RAW FASTQC
########################################

echo "Running FastQC (raw)"

fastqc -t $THREADS \
$FASTQ_DIR/${SAMPLE}_1.fastq.gz \
$FASTQ_DIR/${SAMPLE}_2.fastq.gz \
-o $QC_DIR

multiqc $QC_DIR -o $QC_DIR

########################################
# 2. TRIMMOMATIC
########################################

echo "Running Trimmomatic"

java -jar $TRIMMOMATIC_JAR PE -threads $THREADS \
$FASTQ_DIR/${SAMPLE}_1.fastq.gz \
$FASTQ_DIR/${SAMPLE}_2.fastq.gz \
$TRIM_DIR/${SAMPLE}_1.trimmed.fastq.gz \
$TRIM_DIR/${SAMPLE}_1.unpaired.fastq.gz \
$TRIM_DIR/${SAMPLE}_2.trimmed.fastq.gz \
$TRIM_DIR/${SAMPLE}_2.unpaired.fastq.gz \
ILLUMINACLIP:$ADAPTERS:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
2>&1 | tee $LOG_DIR/${SAMPLE}_trimmomatic.log

########################################
# 3. ALIGNMENT
########################################

echo "Running HISAT2"

hisat2 -p $THREADS \
-x $GENOME_INDEX \
-1 $TRIM_DIR/${SAMPLE}_1.trimmed.fastq.gz \
-2 $TRIM_DIR/${SAMPLE}_2.trimmed.fastq.gz \
--summary-file $LOG_DIR/${SAMPLE}_hisat2_summary.txt \
2> $LOG_DIR/${SAMPLE}_hisat2.log \
| samtools sort -@ $THREADS -o $ALIGN_DIR/${SAMPLE}.sorted.bam

samtools index $ALIGN_DIR/${SAMPLE}.sorted.bam

########################################
# 4. FEATURECOUNTS
########################################

echo "Running featureCounts"

featureCounts \
-T $THREADS \
-p -B -C \
-a $GTF \
-o $COUNT_DIR/${SAMPLE}_counts.txt \
$ALIGN_DIR/${SAMPLE}.sorted.bam \
2>&1 | tee $LOG_DIR/${SAMPLE}_featurecounts.log

########################################
# 5. FLAGSTAT
########################################

echo "Running samtools flagstat"

samtools flagstat $ALIGN_DIR/${SAMPLE}.sorted.bam \
> $LOG_DIR/${SAMPLE}_flagstat.txt

echo "======================================"
echo "Pipeline finished for $SAMPLE"
echo "======================================"
