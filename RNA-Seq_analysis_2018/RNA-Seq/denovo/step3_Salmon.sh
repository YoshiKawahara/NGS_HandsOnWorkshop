#!/bin/bash
ToolDir=$HOME/RNA-Seq/tool
Trinity_bin=$ToolDir/trinityrnaseq-Trinity-v2.8.3
Samtools_bin=$ToolDir/samtools-1.9/bin
Salmon_bin=$ToolDir/salmon-0.11.3-linux_x86_64/bin
R_bin=/work/NGSworkshop/tool/R-3.5.1/bin

export PATH=$Trinity_bin:$Samtools_bin:$Salmon_bin:$R_bin:$PATH

# Step3. Estimating Transcript Abundance by Salmon
$Trinity_bin/util/align_and_estimate_abundance.pl --transcripts Trinity_out/Trinity.fasta \
  --seqType fq --SS_lib_type RF --samples_file sample_info.txt \
  --est_method salmon \
  --thread_count 1 --trinity_mode --prep_reference

# Build Transcript and Gene Expression Matrices
$Trinity_bin/util/abundance_estimates_to_matrix.pl \
  --est_method salmon \
  --gene_trans_map Trinity_out/Trinity.fasta.gene_trans_map \
  --name_sample_by_basedir --out_prefix salmon \
  Day_rep1/quant.sf Day_rep2/quant.sf Day_rep3/quant.sf \
  Night_rep1/quant.sf Night_rep2/quant.sf Night_rep3/quant.sf
