# What is the goal of this project?


## Step 1: Downloading our fastq files.

## Step 2: Getting organized to set ourselves up for success throughout our workflow
```
# Make directories for file organization, change into directory for raw files. Move raw files into the correct directory.

mkdir fastq 
cd fastqfiles
mkdir raw
mkdir logs
mkdir bowtie2
mkdir kraken2
cd raw
```
## Step 3: Running FastQC of your raw data

# FastQC of raw data. This will evaluate the baseline state of our DNA so that we can make decisions about how to clean it later in our workflow.
# Enter interactive mode on a compute node (from where you are)
srun --pty bash

#Load FastQC
module load fastqc

#Confirm that FastQC is available and see options
fastqc -h

# Create a directory for the reports (Once per dataset)
mkdir -p fastqc_out

#Run FastQC. The -o line tells the program to put the .html and .zip outputs into your output directory
fastqc -o fastqc_out SRR6996006.sra_1.fastq.gz SRR6996006.sra_2.fastq.gz
```
The files produced should be:
yourfile_fastqc.html
yourfile_fastqc.zip
```
## Step 4: Run Trimmomatic on the files.
Running trimmomatic is necessary to remove adapters and to trim low quality bases. The quality of reads often tails off at the end of our reads (as evidenced from our initial FastQC analysis. Trimmomatic will remove those while also getting rid of reads that are simply too short for our future needs.

## Step 5: FastQC of clean data:

```
Enter interactive mode on a compute node (from where you are)
srun --pty bash
module load fastqc
fastqc -h
mkdir -p fastqc_out
fastqc -o fastqc_out SRR*
```
Additionally, the trimmed files need to be put in the bucket after they're checked with FastQC.

**Result Files in ReadMe**

## Step 6: Align Cleaned Reads to Human Genome with Bowtie2 and Samtools

# Build human genome index
```
#Download human genome 
wget https://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/hg38/GRCh38.fa.gz

#unzip genome
gunzip GRCh38.fa.gz

#enter compute node to build index 
srun --pty --cpus-per-task=8 --mem=16G --time=02:00:00 bash

#Load module and build index
module load bowtie2
bowtie2-build GRCh38.fa /home/jm3448/Final/bowtie2/GRCh38
```
# Run Bowtie2 and Samtools
```
#Bowtie2 script
#!/bin/bash
#SBATCH --job-name=human_filter
#SBATCH --output=bowtie2_%j.out
#SBATCH --error=bowtie2_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jm3448@georgetown.edu
#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH --mem=16G

module load bowtie2
module load samtools

INDEX="/home/jm3448/Final/bowtie2/GRCh38"
R1="/home/jm3448/Final/trim_out/Control1_paired.fq.gz"
R2="/home/jm3448/Final/trim_out/Control2_paired.fq.gz"

OUTDIR="/home/jm3448/Final/human_filter"
mkdir -p $OUTDIR
cd $OUTDIR

bowtie2 -p 8 -x $INDEX -1 $R1 -2 $R2 | samtools view -bS - > aligned.bam

samtools view -b -f 12 -F 256 aligned.bam > unmapped.bam

samtools fastq \
-1 nonhuman_R1.fastq.gz \
-2 nonhuman_R2.fastq.gz \
unmapped.bam
```
The files produced should be: 
aligned.bam	(all reads--mapped + unmapped)	
unmapped.bam	(non-human reads)
nonhuman_R1.fastq.gz	(filtered forward reads)	
nonhuman_R2.fastq.gz	(filtered reverse reads)	

## Step 7: Classify non-human reads against a reference genome database with Kraken2

# Install Kraken database on login node
```
conda create -n kraken_env -c bioconda kraken2
conda activate kraken_env
```

# Run Kraken2
```
#Kraken2 script
#!/bin/bash
#SBATCH --job-name=kraken2
#SBATCH --output=kraken2_%j.out
#SBATCH --error=kraken2_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jm3448@georgetown.edu
#SBATCH --cpus-per-task=4
#SBATCH --time=08:00:00
#SBATCH --mem=32G

source ~/.bashrc
conda activate kraken_env

DB="/home/jm3448/Final/kraken_db"
R1="/home/jm3448/Final/human_filter/nonhuman_R1.fastq.gz"
R2="/home/jm3448/Final/human_filter/nonhuman_R2.fastq.gz"
OUTDIR="/home/jm3448/Final/kraken_out"

mkdir -p $OUTDIR

kraken2 \
--db $DB \
--threads 4 \
--paired \
--memory-mapping \
$R1 $R2 \
--report $OUTDIR/kraken_report.txt \
--output $OUTDIR/kraken_output.txt
```
The file produced should be: 
kraken_report.txt

# Upload Kraken2 report to the bucket
```
gsutil cp /home/jm3448/Final/kraken_out/kraken_report.txt \
gs://gu-biology-dept-class/jm3448/
```
## Analyze Diversity Metrics of Kraken Report using RStudio

# R Analysis 
```
setwd("~/Downloads")

kraken <- read.table(
  "kraken_report.txt",
  sep="\t",
  fill=TRUE,
  quote="",
  comment.char=""
)

colnames(kraken) <- c("percent","reads","taxReads","rank","taxID","name")
kraken$name <- trimws(kraken$name)

genus <- subset(kraken, rank == "G")

abundance <- genus$reads
names(abundance) <- genus$name
```
# Produce Diversity Metrics
```
library(vegan)

shannon <- diversity(abundance, index="shannon")
simpson <- diversity(abundance, index="simpson")
richness <- specnumber(abundance)
```
# Visualize Data
```
rel_abundance <- abundance / sum(abundance)
top10 <- sort(rel_abundance, decreasing=TRUE)[1:10]

par(mar=c(10,4,4,2))

barplot(
  top10,
  las=2,
  col="darkgreen",
  main="Top 10 Genera (Relative Abundance)",
  ylab="Proportion"
)
```
**Result Files in ReadMe**

