#!/bin/bash

DataDir=$HOME/RNA-Seq/data
ToolDir=$HOME/RNA-Seq/tool
STAR_bin=$ToolDir/STAR-2.6.1b/bin/Linux_x86_64

export PATH=$STAR_bin:$PATH

### Step2. Alignment of RNA-Seq reads by STAR
GenomeDir=STAR_index
Read1_list=$DataDir/rice_D_rep1_r1.org.fastq.gz,$DataDir/rice_D_rep2_r1.org.fastq.gz,$DataDir/rice_D_rep3_r1.org.fastq.gz,$DataDir/rice_N_rep1_r1.org.fastq.gz,$DataDir/rice_N_rep2_r1.org.fastq.gz,$DataDir/rice_D_rep3_r1.org.fastq.gz
Read2_list=$DataDir/rice_D_rep1_r2.org.fastq.gz,$DataDir/rice_D_rep2_r2.org.fastq.gz,$DataDir/rice_D_rep3_r2.org.fastq.gz,$DataDir/rice_N_rep1_r2.org.fastq.gz,$DataDir/rice_N_rep2_r2.org.fastq.gz,$DataDir/rice_D_rep3_r2.org.fastq.gz

STAR --genomeDir $GenomeDir --readFilesCommand "gunzip -c" \
  --readFilesIn $Read1_list $Read2_list \
  --alignIntronMin 20 --alignIntronMax 10000 \
  --outSAMtype BAM SortedByCoordinate \
  --runThreadN 1
