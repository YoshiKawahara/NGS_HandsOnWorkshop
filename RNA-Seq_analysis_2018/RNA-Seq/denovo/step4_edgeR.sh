#!/bin/bash
ToolDir=$HOME/RNA-Seq/tool
Trinity_bin=$ToolDir/trinityrnaseq-Trinity-v2.8.3
R_bin=/work/NGSworkshop/tool/R-3.5.1/bin

export PATH=$Trinity_bin:$R_bin:$PATH

# Step 4. Differential Expression Analysis by edgeR
$Trinity_bin/Analysis/DifferentialExpression/run_DE_analysis.pl \
  --matrix salmon.gene.counts.matrix \
  --method edgeR \
  --samples_file sample_info.txt \
  --output edgeR_out
