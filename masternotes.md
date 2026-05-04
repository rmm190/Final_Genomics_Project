# What is the goal of this project?
The goal of this project is to perform a bioinformatics workflow, similar to a study performed by Li et al. 2021, in order to analyze how the diversity and compostion of the gut microbiome varies with COVID-19 disease severity. To accomplish this, 3 representative fecal samples were studied: one healthy control, and two COVID-19 patients of mild and severe disease severity.  We aim to replicate and adjust the workflow of Li et al. 2021 in order to both evaluate these patterns in gut microbial composition associated with COVID 19 infection and to become more comfortable using common bioinformatics tools. 

**All of the following code is representative. The same protocol and workflow was applied to all of the samples tested.***

## Step 1: Downloading our fastq files.
#We placed all fastq files obtained from the Li et al. study in the bucket.
```
gsutil cp gs://gu-biology-dept-class/jm3448/Control1.fastq.gz
```
## Step 2: Getting organized to set ourselves up for success throughout our workflow
# Make directories for file organization, change into directory for raw files. Move raw files into the correct directory.

```
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
```
# Enter interactive mode on a compute node (from where you are)
srun --pty bash

#Load FastQC
module load fastqc

#Confirm that FastQC is available and see options
fastqc -h

# Create a directory for the reports (Once per dataset)
mkdir -p fastqc_out

#Run FastQC. The -o line tells the program to put the .html and .zip outputs into your output directory
fastqc -o fastqc_out control.sra_1.fastq.gz control.sra_2.fastq.gz
```
The files produced should be:
yourfile_fastqc.html
yourfile_fastqc.zip

## Step 4: Run Trimmomatic on the files.
Running trimmomatic is necessary to remove adapters and to trim low quality bases. The quality of reads often tails off at the end of our reads (as evidenced from our initial FastQC analysis). Trimmomatic will remove those while also getting rid of reads that are simply too short for our future needs.
```
#SBATCH --mail-type=END,FAIL --mail-user=rmm190@georgetown.edu
#SBATCH --job-name="trim_S1"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=03:00:00
#SBATCH --mem=10G

#Load Trimmomatic module ("aliases" needed for GU HPC setup here)
shopt -s expand_aliases
module load trimmomatic

#Define paths and variables
R1=/home/rmm190/final_project/raw_reads/Mild1.fastq.gz
R2=/home/rmm190/final_project/raw_reads/Mild2.fastq.gz
OUTDIR=/home/rmm190/final_project/trim_out
mkdir -p $OUTDIR

trimmomatic PE -threads $SLURM_CPUS_PER_TASK \
$R1 $R2 \
$OUTDIR/Mild_1_paired.fq.gz $OUTDIR/Mild_1_unpaired.fq.gz \
$OUTDIR/Mild_2_paired.fq.gz $OUTDIR/Mild_2_unpaired.fq.gz \
ILLUMINACLIP:$HOME/rmm190/TruSeq3-PE.fa:2:30:10 \
LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:75

```
## Step 5: FastQC of clean data:

```
Enter interactive mode on a compute node (from where you are)
srun --pty bash
module load fastqc
fastqc -h
mkdir -p fastqc_out
fastqc -o fastqc_out Control*
```
Additionally, the trimmed files need to be put in the bucket after they're checked with FastQC.


## Step 6: Align Cleaned Reads to Human Genome with Bowtie2 and Samtools
The purpose of this step is to match human DNA present within the samples to a humman reference genome to remove it from our analyses. We're focused on the gut microbiome, and having human DNA present will inhibit our ability to analyze exclusively the gut microbiome. Bowtie2 is specifically for aligning the human reads to the human reference genome, and Samtools is for removing those human reads. At the end of this step, we want to keep the unmapped reads which represent the microbial DNA. 

# Build human genome index
This step produces a searchable genome that enables a more efficinet search of the human genome. 
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
The purpose of this step is to identify the organisms that are present in our sample. Kraken2 works specifically by breaking reads into k-mers and then comparing these pieces against a known reference database to match reads to taxonomies present.

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
R studio is an important tool for visualizing different diversity metrics within our cleaned reads. This will help us obtain a better understanding of the gut microbiome composition in both healthy patients and COVID 19 patients across the two tested levels of severity: Mild and Severe. Metrics analyzed in R include relative abundance, richness, Shannon diversity index and Simpson index. These metrics are slightly different from those analyzed in the original Li et al., study and can hopefully provide a more nuanced perspective of gut microbiome composition following COVID-19 infection of different severity.  

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

# Install package (run once if needed)
install.packages("tidyverse")

# Load library
library(tidyverse)

# ---------------------------
# Step 1: Create data
# ---------------------------

control <- data.frame(
  Genus = c("Bacteroides","Phocaeicola","Faecalibacterium","Agathobacter",
            "Alistipes","Parabacteroides","Klebsiella","Dialister",
            "Mediterraneibacter","Roseburia"),
  Abundance = c(0.30,0.19,0.14,0.048,0.042,0.036,0.030,0.025,0.022,0.018),
  Group = "Control"
)

severe <- data.frame(
  Genus = c("Bacteroides","Phocaeicola","Roseburia","Phascolarctobacterium",
            "Parabacteroides","Faecalibacterium","Sinaaca","Lachnospira",
            "Bifidobacterium","Haemophilus"),
  Abundance = c(0.57,0.29,0.025,0.025,0.015,0.015,0.010,0.008,0.007,0.005),
  Group = "Severe"
)

mild <- data.frame(
  Genus = c("Segatella","Bacteroides","Phocaeicola","Faecalibacterium",
            "Parabacteroides","Alistipes","Dialister","Roseburia",
            "Agathobacter","Butyricimonas"),
  Abundance = c(0.81,0.10,0.025,0.020,0.015,0.012,0.010,0.009,0.008,0.007),
  Group = "Mild"
)

# ---------------------------
# Step 2: Combine datasets
# ---------------------------

df <- bind_rows(control, severe, mild)

# Fill missing genera with 0
df_complete <- df %>%
  complete(Genus, Group, fill = list(Abundance = 0))

# ---------------------------
# Step 3: Order by Control abundance
# ---------------------------

genus_order <- df_complete %>%
  filter(Group == "Control") %>%
  arrange(desc(Abundance)) %>%
  pull(Genus)

# Add genera not present in Control to the end
genus_order <- c(genus_order,
                 setdiff(unique(df_complete$Genus), genus_order))

# Apply ordering
df_complete$Genus <- factor(df_complete$Genus, levels = genus_order)

# ---------------------------
# Step 4: Plot
# ---------------------------

ggplot(df_complete, aes(x = Genus, y = Abundance, fill = Group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Control" = "green",
                               "Severe" = "red",
                               "Mild" = "blue")) +
  labs(title = "Stacked Relative Abundance of Top Genera",
       x = "Genus",
       y = "Proportion") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
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
# Install if needed
install.packages("tidyverse")

# Load
library(tidyverse)

# ---------------------------
# Create data
# ---------------------------

df <- data.frame(
  Metric = rep(c("Richness", "Shannon Diversity", "Simpson Index"), each = 3),
  Group = rep(c("Control", "Mild", "Severe"), times = 3),
  Value = c(
    2662, 56, 1507,
    2.627657, 0.9137039, 1.27922,
    0.8408519, 0.782076, 0.5761584
  )
)

# ---------------------------
# Plot
# ---------------------------

ggplot(df, aes(x = Group, y = Value, fill = Group)) +
  geom_bar(stat = "identity", width = 0.6) +
  
  facet_wrap(~Metric, scales = "free_y") +
  
  scale_fill_manual(values = c(
    "Control" = "darkgreen",
    "Mild" = "blue",
    "Severe" = "red"
  )) +
  
  labs(
    title = "Microbial Diversity Metrics Across Groups",
    x = "",
    y = "Value"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )
```
**Result Files in ReadMe**

