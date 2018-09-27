#!/bin/bash
ToolDir=$HOME/RNA-Seq/tool
Trinity_bin=$ToolDir/trinityrnaseq-Trinity-v2.8.3
Samtools_bin=$ToolDir/samtools-1.9/bin
Jellyfish_bin=$ToolDir/jellyfish-2.2.10/bin
Salmon_bin=$ToolDir/salmon-0.11.3-linux_x86_64/bin
Bowtie2_bin=$ToolDir/bowtie2-2.3.4.2-linux-x86_64

export PATH=$Trinity_bin:$Samtools_bin:$Jellyfish_bin:$Salmon_bin:$Bowtie2_bin:$PATH

# Step2. De novo transcriptome assemble by Trinity
Trinity --seqType fq --max_memory 4G --CPU 1 --SS_lib_type RF --samples_file sample_info.txt --output Trinity_out
