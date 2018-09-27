#!/bin/bash

DataDir=$HOME/RNA-Seq/data
ToolDir=$HOME/RNA-Seq/tool
HISAT2_bin=$ToolDir/hisat2-2.1.0
Samtools_bin=$ToolDir/samtools-1.9

export PATH=$HISAT2_bin:$Samtools_bin:$PATH

### Step3. Align reads to the reference genome by HISAT2
HISAT2_COMMON_PARAM="--threads 1 --min-intronlen 20 --max-intronlen 10000 --dta --rna-strandness RF -x genome"

for DATASET in rice_D_rep1 rice_D_rep2 rice_D_rep3 rice_N_rep1 rice_N_rep2 rice_N_rep3
do
  hisat2 $HISAT2_COMMON_PARAM -1 ${DATASET}_r1.pe.fastq.gz -2 ${DATASET}_r2.pe.fastq.gz -S ${DATASET}.sam
  samtools sort -o ${DATASET}.bam ${DATASET}.sam
  samtools index ${DATASET}.bam
  rm ${DATASET}.sam
done
