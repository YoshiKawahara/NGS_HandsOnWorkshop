#!/bin/bash

mkdir tool
cd tool
ToolDir=`pwd`

echo "START:" `date`
### Tools for "useref" analysis
echo "installing FastQC..."
wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip
unzip fastqc_v0.11.7.zip
chmod +x FastQC/fastqc
rm fastqc_v0.11.7.zip
echo "installing Trimmomatic..."
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip
unzip Trimmomatic-0.38.zip
rm Trimmomatic-0.38.zip
echo "installing HISAT2..."
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip
unzip hisat2-2.1.0-Linux_x86_64.zip
rm hisat2-2.1.0-Linux_x86_64.zip
echo "installing StringTie..."
wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.4d.Linux_x86_64.tar.gz
tar xfz stringtie-1.3.4d.Linux_x86_64.tar.gz
rm stringtie-1.3.4d.Linux_x86_64.tar.gz
echo "installing Samtools..."
wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
tar xfj samtools-1.9.tar.bz2
rm samtools-1.9.tar.bz2
cd samtools-1.9
make
make prefix=`pwd` install
cd ..

### Tools for "denovo" analysis
echo "installing CMake..."
wget https://cmake.org/files/v3.12/cmake-3.12.1-Linux-x86_64.tar.gz
tar xfz cmake-3.12.1-Linux-x86_64.tar.gz
rm cmake-3.12.1-Linux-x86_64.tar.gz
PATH=$ToolDir/cmake-3.12.1-Linux-x86_64/bin:$PATH
export PATH
echo "installing Trinity..."
wget https://github.com/trinityrnaseq/trinityrnaseq/archive/Trinity-v2.8.3.tar.gz
tar xfz Trinity-v2.8.3.tar.gz
rm Trinity-v2.8.3.tar.gz
cd trinityrnaseq-Trinity-v2.8.3
make
cd ..
echo "installing Bwotie2..."
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.4.2/bowtie2-2.3.4.2-linux-x86_64.zip
unzip bowtie2-2.3.4.2-linux-x86_64.zip
rm bowtie2-2.3.4.2-linux-x86_64.zip
echo "installing Jellyfish..."
wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.10/jellyfish-2.2.10.tar.gz
tar xfz jellyfish-2.2.10.tar.gz
rm jellyfish-2.2.10.tar.gz
cd jellyfish-2.2.10
./configure --prefix=`pwd`
make
make install
cd ..
echo "installing Salmon..."
wget https://github.com/COMBINE-lab/salmon/releases/download/v0.11.3/salmon-0.11.3-linux_x86_64.tar.gz
tar xfz salmon-0.11.3-linux_x86_64.tar.gz
rm salmon-0.11.3-linux_x86_64.tar.gz
echo "installing BLAST+..."
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST//ncbi-blast-2.7.1+-x64-linux.tar.gz
tar xfz ncbi-blast-2.7.1+-x64-linux.tar.gz
rm ncbi-blast-2.7.1+-x64-linux.tar.gz

### Tools for "varcall" analysis
echo "installing STAR..."
wget https://github.com/alexdobin/STAR/archive/2.6.1b.tar.gz
tar xfz 2.6.1b.tar.gz
rm 2.6.1b.tar.gz
echo "installing Picard..."
mkdir picard-2.18.12
cd picard-2.18.12
wget https://github.com/broadinstitute/picard/releases/download/2.18.12/picard.jar
cd ..
echo "installing GATK..."
wget https://github.com/broadinstitute/gatk/releases/download/4.0.8.1/gatk-4.0.8.1.zip
unzip gatk-4.0.8.1.zip
rm gatk-4.0.8.1.zip

echo "END:" `date`
