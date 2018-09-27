#!/bin/bash

# Step1. Add suffix (/1 or /2) to original FASTQ header for Trinity
for DATASET in rice_D_rep1 rice_D_rep2 rice_D_rep3 rice_N_rep1 rice_N_rep2 rice_N_rep3
do
    echo "converting ${DATASET}_r1 ..."
    zcat ../data/${DATASET}_r1.org.fastq.gz | awk '{ if (NR%4==1) { print $1"/1" } else { print } }' | gzip -c > ${DATASET}_r1.trinity.fastq.gz
    echo "converting ${DATASET}_r2 ..."
    zcat ../data/${DATASET}_r2.org.fastq.gz | awk '{ if (NR%4==1) { print $1"/2" } else { print } }' | gzip -c > ${DATASET}_r2.trinity.fastq.gz
done
