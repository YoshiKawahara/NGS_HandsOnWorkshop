#!/bin/bash

DataDir=$HOME/RNA-Seq/data
ToolDir=$HOME/RNA-Seq/tool
Picard_bin=$ToolDir/picard-2.18.12
Samtools_bin=$ToolDir/samtools-1.9
GATK_bin=$ToolDir/gatk-4.0.8.1

export PATH=$Samtools_bin:$GATK_bin:$PATH

### Step3. Processing of alignment (BAM) file
java -jar $Picard_bin/picard.jar AddOrReplaceReadGroups \
  I=Aligned.sortedByCoord.out.bam \
  O=Aligned.sortedByCoord.RG.out.bam \
  SO=coordinate \
  RGID=Koshihikari RGLB=TruSeq_RNA_stranded RGPL=illumina RGPU=HiSeq2000 RGSM=Koshihikari

java -jar $Picard_bin/picard.jar MarkDuplicates \
  I=Aligned.sortedByCoord.RG.out.bam \
  O=Aligned.sortedByCoord.RG.MD.out.bam \
  CREATE_INDEX=true \
  VALIDATION_STRINGENCY=SILENT \
  M=Aligned.sortedByCoord.RG.MD.metrics

samtools faidx $DataDir/genome.fa

java -jar $Picard_bin/picard.jar CreateSequenceDictionary \
  R=$DataDir/genome.fa \
  O=$DataDir/genome.dict

gatk SplitNCigarReads \
   -R $DataDir/genome.fa \
   -I Aligned.sortedByCoord.RG.MD.out.bam \
   -O Aligned.sortedByCoord.RG.MD.SplitN.out.bam
