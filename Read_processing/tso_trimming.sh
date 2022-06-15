#!/bin/bash
#SBATCH -N 1 -n 1 -c 2
#SBATCH --mem=20000
#SBATCH -t 18:00:00
#SBATCH -J tso_trimming
#SBATCH --mail-type=ALL
#SBATCH -e job-%J.err
#SBATCH -o job-%J.out

echo "Execuuting file" $0
echo "sample R2" $1
echo "Experiment" $2
echo "subarray" $3
#echo "sample R1" $4

#File name
sn=$1
fn=$2
qn=$3

#Go to the folder where fastq files are
cd Arabidopsis/10_45_45_all_samples/reads/Read2/
#Unzip the files
gunzip raw/${sn}_R2.fastq.gz

#TSO trimming
cutadapt --cores=0 \
    -g ""XAAGCAGTGGTATCAACGCAGAGTACATGGG";max_error_rate="0.1"" \
    --overlap 5 \
    -o Arabidopsis/10_45_45_all_samples/reads/Read2/tso_trimmed/${fn}_${qn}_trimmed_R2.fastq \
    -n 2 raw/${sn}_R2.fastq

echo "TSO trimming done!"

#Extract headers from read 2, read longer than 0 bp
Softwares/seqtk/seqtk comp Arabidopsis/10_45_45_all_samples/reads/Read2/tso_trimmed/${fn}_${qn}_trimmed_R2.fastq | awk '{ if ($2 > 0) { print} }' | cut --fields 1 > Arabidopsis/10_45_45_all_samples/reads/Read2/tso_trimmed/${sn}_headers.txt

echo "Headers extracted from R2!"

#Filter R2
Softwares/seqtk/seqtk subseq Arabidopsis/10_45_45_all_samples/reads/Read2/tso_trimmed/${fn}_${qn}_trimmed_R2.fastq Arabidopsis/10_45_45_all_samples/reads/Read2/tso_trimmed/${sn}_headers.txt > Arabidopsis/10_45_45_all_samples/reads/Read2/tso_trimmed/${fn}_${qn}_trimmed_filtered_R2.fastq

echo "Read 2 processed!"

#move to folder with R1
cd Arabidopsis/10_45_45_all_samples/reads/Read1

#Unzip R1
gunzip *.gz

#Filter R1
Softwares/seqtk/seqtk subseq ${sn}_R1.fastq Arabidopsis/10_45_45_all_samples/reads/Read2/tso_trimmed/${sn}_headers.txt > Arabidopsis/10_45_45_all_samples/reads/Read1/filtered/${fn}_${qn}_filtered_R1.fastq

echo "Read1 filtered!"

#zip the original R1
pigz ${sn}_R1.fastq

echo "Read1 done!!"
