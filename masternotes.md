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
