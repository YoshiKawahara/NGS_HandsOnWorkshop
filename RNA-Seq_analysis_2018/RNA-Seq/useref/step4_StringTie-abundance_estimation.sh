#!/bin/bash

DataDir=$HOME/RNA-Seq/data
ToolDir=$HOME/RNA-Seq/tool
StringTie_bin=$ToolDir/stringtie-1.3.4d.Linux_x86_64

export PATH=$StringTie_bin:$PATH

### Step4. Estimate transcript abundances
STRINGTIE_COMMON_PARAM="-e -B"
for DATASET in rice_D_rep1 rice_D_rep2 rice_D_rep3 rice_N_rep1 rice_N_rep2 rice_N_rep3 
do
  stringtie $STRINGTIE_COMMON_PARAM -G $DataDir/annotation.gtf -o ballgown/${DATASET}/${DATASET}.gtf ${DATASET}.bam
done

