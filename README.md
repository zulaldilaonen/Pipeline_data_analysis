# RNA-Seq Analysis Pipeline (Single Sample)

This repository contains a reproducible RNA-Seq analysis workflow for a paired-end human sample (example: SRR14513977).

The pipeline performs:

1. Raw FASTQ quality control (FastQC)
2. Adapter and quality trimming (Trimmomatic)
3. Post-trimming QC (FastQC + MultiQC)
4. Alignment to GRCh38 reference genome (HISAT2)
5. BAM sorting and indexing (samtools)
6. Gene-level read counting (featureCounts)

Reference Genome: GRCh38 (Primary Assembly)  
Annotation: GENCODE v45  

---

# 1. System Requirements

Recommended environment:
- Ubuntu 20.04+ or WSL2
- Minimum 16 GB RAM
- 50–100 GB free disk space

Update system:

```bash
sudo apt update && sudo apt upgrade -y
```
---

# 2. Install Required Tools

Install all required bioinformatics tools:

```bash 
sudo apt install -y \
fastqc \
multiqc \
hisat2 \
samtools \
subread \
trimmomatic \
wget \
gzip \
default-jre
```
Explanation of tools:

FastQC → Quality control

MultiQC → QC summary report

Trimmomatic → Adapter and quality trimming

HISAT2 → Alignment

samtools → BAM processing

featureCounts (from subread) → Gene counting

Verify installation:

```bash 
fastqc --version
hisat2 --version
samtools --version
featureCounts -v

```
---

# 3. Download Reference Genome

Create reference directory:

```bash 
mkdir reference
cd reference

```
Download GRCh38 genome and GENCODE v45 annotation:
```bash 
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/GRCh38.primary_assembly.genome.fa.gz

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz

```

Uncompress:

```bash
gunzip *.gz
```
Build HISAT2 index:
```bash
hisat2-build GRCh38.primary_assembly.genome.fa GRCh38_index
```

This will generate files:
```bash
ls
```
Output like this:
GRCh38_index.1.ht2
GRCh38_index.2.ht2
...
GRCh38_index.8.ht2

# 4. Create Adapter File
Create adapters.fa inside reference directory:
```bash
vim adapters.fa
```
Paste:
```bash
>Adapter1
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
>Adapter2
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

```
Save and exit.(esc, then write :wq)


# 5. Download FASTQ Data
Return to main directory:
```bash
cd ..
mkdir fastq
cd fastq

```
Download paired-end FASTQ from ENA:
```bash
wget https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/077/SRR14513977/SRR14513977_1.fastq.gz

wget https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR145/077/SRR14513977/SRR14513977_2.fastq.gz

```
Return to main directory:
```bash
cd ..

```

# 6. Run the Pipeline
Make script executable:

```bash
chmod +x pipeline.sh

```
Run analysis:
```bash
./pipeline.sh SRR14513977

```
# 7. Output Structure

After completion, the following folders will be created automatically:

fastq/          → Raw FASTQ files
trimmed/        → Trimmed reads
alignment/      → Sorted BAM + index
counts/         → Gene count table
QC/             → FastQC & MultiQC reports
logs/           → Tool log files
reference/      → Genome and annotation files


# 8. Pipeline Steps Explained
## Step 1 – Raw Quality Control

FastQC analyzes:

-Per-base quality

-GC content

-Adapter contamination

-Sequence duplication

## Step 2 – Trimming

Trimmomatic removes:

-Illumina adapters

-Low-quality bases

-Reads shorter than 36 bp

## Step 3 – Alignment

-HISAT2 aligns cleaned reads to GRCh38 reference genome.

##Step 4 – Gene Counting

-featureCounts assigns aligned reads to genes using GENCODE v45 annotation.

# 9. Notes

-Large files (FASTQ, BAM, genome index) are not tracked in GitHub.

-The pipeline assumes paired-end RNA-Seq data.

-For strand-specific libraries, featureCounts parameters may need adjustment.
# 10. Example Expected Metrics (Typical RNA-Seq)

-Trimmed paired-end reads: >70%

-Alignment rate: 70–90%

-Gene assignment rate: 60–80%

Lower values may indicate:

-Library type mismatch

-rRNA contamination

-Strand mis-specification

-Low input quality
