#!/bin/bash
#SBATCH -N 1 -n 1 -c 3
#SBATCH --mem=50000
#SBATCH -t 10:00:00
#SBATCH -J EXP_mm
#SBATCH --mail-type=ALL
#SBATCH -e job-%J.err
#SBATCH -o job-%J.out

module load samtools/1.6

echo "Execuuting file" $0
echo "Sample ID" $1

#Folder name
sn=$1

#Create folders for run where R2 is multimapped and discarded reads are processed.
mkdir temp/${sn}_bacteria
mkdir data/${sn}_bacteria

# Reads and trimming
FW=Arabidopsis/10_45_45_all_samples/reads/Read1/spatial/${sn}_filtered_R11.fastq
RV=Arabidopsis/10_45_45_all_samples/reads/Read2/tso_trimmed/${sn}_trimmed_filtered_R2.fastq

# References for mapping, annotation and ribo-filtering
MAP=Genomes_transcriptomes_indexies/Arabidopsis/Arabidopsis_genome_mCitrine/StarIndex_Arabidopsis_Pseudomonas
ANN=Genomes_transcriptomes_indexies/Arabidopsis/Arabidopsis_genome_mCitrine/TAIR10_genes.gtf

# Barcodes settings
ID=Genomes_transcriptomes_indexies/omni-v1_coordinates.txt

# Output folder and experiment name
OUTPUT=Arabidopsis/10_45_45_all_samples/data/${sn}_bacteria
EXP=${sn}_bacteria
TEMP=Arabidopsis/10_45_45_all_samples/temp/${sn}_bacteria

#ST pipeline to allow Arabidopsis multimapping
st_pipeline_run.py \
  --ids $ID \
  --ref-map $MAP \
  --expName $EXP \
  --ref-annotation $ANN \
  --mapping-threads 25 \
  --allowed-missed 1 \
  --allowed-kmer 4 \
  --umi-allowed-mismatches 2 \
  --umi-start-position 16 \
  --umi-end-position 28 \
  --output-folder $OUTPUT \
  --verbose \
  --remove-polyA 15 \
  --remove-polyT 15 \
  --remove-polyC 15 \
  --remove-polyG 15 \
  --compute-saturation \
  --keep-discarded-files \
  --disable-annotation \
  --temp-folder $TEMP \
  --no-clean-up \
  --log-file $OUTPUT/${EXP}_log.txt \
  $FW $RV

#Rename mapping_discarded to mapped --> bacterial reads
mv $TEMP/mapping_discarded.bam $TEMP/mapped.bam

#ST pipeline to demultiplex bacterial reads
st_pipeline_run.py \
  --ids $ID \
  --ref-map $MAP \
  --expName $EXP \
  --ref-annotation $ANN \
  --mapping-threads 15 \
  --allowed-missed 1 \
  --allowed-kmer 4 \
  --umi-allowed-mismatches 2 \
  --umi-start-position 16 \
  --umi-end-position 28 \
  --output-folder $OUTPUT \
  --verbose \
  --remove-polyA 15 \
  --remove-polyT 15 \
  --remove-polyC 15 \
  --remove-polyG 15 \
  --disable-multimap \
  --compute-saturation \
  --keep-discarded-files \
  --disable-trimming \
  --disable-umi \
  --disable-mapping \
  --disable-annotation \
  --temp-folder $TEMP \
  --no-clean-up \
  --log-file $OUTPUT/${EXP}_log.txt \
  $FW $RV

#Clip the header information to a new file
samtools view $TEMP/demultiplexed_matched.bam | cut -f 1,18,21,22 > $TEMP/${sn}_demultiplex_matched.tsv

#Remove extra from the names
sed -e s/B1:Z://g -e s/B3:Z://g -e s/$'\t'B2:Z:/x/g $TEMP/${sn}_demultiplex_matched.tsv > $TEMP/${sn}_demultiplex_matched2.tsv

#Rename the file as origina,
mv $TEMP/${sn}_demultiplex_matched2.tsv $TEMP/${sn}_demultiplex_matched.tsv

#Collect the reads to a fastq file
samtools bam2fq $TEMP/demultiplexed_matched.bam > $OUTPUT/${sn}_bacteria_R2.fastq

#zip the fastq and tsv files
pigz $OUTPUT/${sn}_bacteria_R2.fastq
pigz $TEMP/${sn}_demultiplex_matched.tsv

#Change conda environment
conda deactivate
conda activate tagGD_demultiplex

#Extract the unmapped headers
samtools view $TEMP/mapped.bam | cut -f 1 > $TEMP/${sn}_unmapped_headers.txt

#Filter the probe sequences based on unmapped headers
Softwares/seqtk/seqtk subseq Arabidopsis/10_45_45_all_samples/reads/Read1/probe/${sn}_filtered_R12.fastq $TEMP/${sn}_unmapped_headers.txt > Arabidopsis/10_45_45_all_samples/reads/Read1/probe/${sn}_discarded_filtered_R12.fastq

#TagGD for probe demultiplexing
taggd_demultiplex.py \
        --k 5 \
        --start-position 2 \
        --max-edit-distance 3 \
        --seed 8946 \
        --no-offset-speedup \
        --multiple-hits-keep-one \
        Arabidopsis/degenerative_probes.tsv \
        Arabidopsis/10_45_45_all_samples/reads/Read1/probe/${sn}_discarded_filtered_R12.fastq Arabidopsis/10_45_45_all_samples/reads/Read1/probe/${sn}

#Extract headers
sed -n '1~4p' Arabidopsis/10_45_45_all_samples/reads/Read1/probe/${sn}_matched.fastq | cut -d \  -f 1,3,4> Arabidopsis/10_45_45_all_samples/reads/Read1/probe/${sn}_probes.txt

#Remove extra parts from the names
sed -e s/B0:Z://g -e s/B1:Z://g Arabidopsis/10_45_45_all_samples/reads/Read1/probe/${sn}_probes.txt > Arabidopsis/10_45_45_all_samples/reads/Read1/probe/${sn}_probe_information.txt

#zip the probe file
pigz Arabidopsis/10_45_45_all_samples/reads/Read1/probe/${sn}_probe_information.txt
