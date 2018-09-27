#!/bin/bash

DataDir=$HOME/RNA-Seq/data
ToolDir=$HOME/RNA-Seq/tool
GATK_bin=$ToolDir/gatk-4.0.8.1

export PATH=$GATK_bin:$PATH

### Step4. Detection of variations by GATK HaplotypeCaller
gatk HaplotypeCaller -R $DataDir/genome.fa -I Aligned.sortedByCoord.RG.MD.SplitN.out.bam \
  -stand-call-conf 20 --dont-use-soft-clipped-bases \
  -O variation.vcf.gz

gatk VariantFiltration -R $DataDir/genome.fa -V variation.vcf.gz \
  -window 20 -cluster 3 \
  --filter-name "QD" --filter "QD < 2.0" \
  --filter-name "FS" --filter "FS > 60.0" \
  -O variation.filter.vcf.gz
