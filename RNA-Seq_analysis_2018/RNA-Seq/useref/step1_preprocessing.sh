#!/bin/bash

DataDir=$HOME/RNA-Seq/data
ToolDir=$HOME/RNA-Seq/tool
FastQC_bin=$ToolDir/FastQC
Trimmomatic_bin=$ToolDir/Trimmomatic-0.38

export PATH=$FastQC_bin:$PATH

### Step1. Preprocessing by FastQC and Trimmomatic
FASTQC_OUTDIR_BEFORE=FastQC_before_preprocess
FASTQC_OUTDIR_AFTER=FastQC_after_preprocess
mkdir $FASTQC_OUTDIR_BEFORE $FASTQC_OUTDIR_AFTER

for DATASET in rice_D_rep1 rice_D_rep2 rice_D_rep3 rice_N_rep1 rice_N_rep2 rice_N_rep3 
do
  # FastQC before Trimmomatic
  fastqc --threads 1 --nogroup  --outdir $FASTQC_OUTDIR_BEFORE --format fastq $DataDir/${DATASET}_r*.org.fastq.gz 
  # Trimmomatic
  java -Xmx4G -Xms2G -jar $Trimmomatic_bin/trimmomatic-0.38.jar PE \
  -threads 1 -phred33 -trimlog Trimmomatic_${DATASET}.log  \
  $DataDir/${DATASET}_r1.org.fastq.gz  $DataDir/${DATASET}_r2.org.fastq.gz \
  ${DATASET}_r1.pe.fastq.gz ${DATASET}_r1.unpe.fastq.gz \
  ${DATASET}_r2.pe.fastq.gz ${DATASET}_r2.unpe.fastq.gz \
  ILLUMINACLIP:$Trimmomatic_bin/adapters/TruSeq3-PE-2.fa:2:30:10 \
  LEADING:15 TRAILING:15 SLIDINGWINDOW:10:15 MINLEN:50
  # FastQC after Trimmomatic
  fastqc --threads 1 --nogroup --outdir $FASTQC_OUTDIR_AFTER --format fastq ${DATASET}_r*.pe.fastq.gz
done
