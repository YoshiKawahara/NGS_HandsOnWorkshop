#!/bin/bash
ToolDir=$HOME/RNA-Seq/tool
Bowtie2_bin=$ToolDir/bowtie2-2.3.4.2-linux-x86_64
Samtools_bin=$ToolDir/samtools-1.9/bin
BLAST_bin=$ToolDir/ncbi-blast-2.7.1+/bin
Trinity_bin=$ToolDir/trinityrnaseq-Trinity-v2.8.3

export PATH=$Bowtie2_bin:$Samtools_bin:$BLAST_bin:$PATH

### Step5. Evalution of de novo transcriptome assemblies by Trinity
# check read representation
bowtie2-build Trinity_out/Trinity.fasta Trinity.fasta
bowtie2 -q --no-unal -k 20 -x Trinity.fasta \
  -1 rice_D_rep1_r1.trinity.fastq.gz,rice_D_rep2_r1.trinity.fastq.gz,rice_D_rep3_r1.trinity.fastq.gz,rice_N_rep1_r1.trinity.fastq.gz,rice_N_rep2_r1.trinity.fastq.gz,rice_N_rep3_r1.trinity.fastq.gz \
  -2 rice_D_rep1_r2.trinity.fastq.gz,rice_D_rep2_r2.trinity.fastq.gz,rice_D_rep3_r2.trinity.fastq.gz,rice_N_rep1_r2.trinity.fastq.gz,rice_N_rep2_r2.trinity.fastq.gz,rice_N_rep3_r2.trinity.fastq.gz \
  | samtools view -Sb -o bowtie2_read_to_contig.bam 

# homology search against rice protein sequences
ln -s ../data/rice_protein.fa
$BLAST_bin/makeblastdb -dbtype prot -in rice_protein.fa
$BLAST_bin/blastx -query Trinity_out/Trinity.fasta -db rice_protein.fa \
  -task blastx-fast -evalue 1e-20 -max_target_seqs 1 -outfmt 6 \
  -out blastx_trinity_to_rice_protein.outfmt6

perl $Trinity_bin/util/analyze_blastPlus_topHit_coverage.pl blastx_trinity_to_rice_protein.outfmt6 Trinity_out/Trinity.fasta rice_protein.fa > blastx_TopHitCoverage.txt
