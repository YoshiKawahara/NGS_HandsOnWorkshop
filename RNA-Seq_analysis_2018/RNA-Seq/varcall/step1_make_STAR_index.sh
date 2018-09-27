#!/bin/bash

DataDir=$HOME/RNA-Seq/data
ToolDir=$HOME/RNA-Seq/tool
STAR_bin=$ToolDir/STAR-2.6.1b/bin/Linux_x86_64

export PATH=$STAR_bin:$PATH

### Step1. Make genome index for STAR
GenomeDir=STAR_index
mkdir $GenomeDir

STAR --runMode genomeGenerate \
  --genomeFastaFiles $DataDir/genome.fa \
  --genomeDir $GenomeDir \
  --sjdbGTFfile $DataDir/annotation.gtf \
  --sjdbOverhang 100 --runThreadN 1
