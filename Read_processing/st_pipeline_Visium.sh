#!/bin/bash
#SBATCH -N 1 -n 1 -c 4
#SBATCH --mem=100000
#SBATCH -t 20:00:00
#SBATCH -J EXP13
#SBATCH --mail-type=ALL
#SBATCH -e job-%J.err
#SBATCH -o job-%J.out

module load samtools/1.6

echo "Execuuting file" $0
echo "Sample ID" $1

#Folder name
sn=$1

mkdir data/${sn}
mkdir temp/${sn}

# Reads
FW=Arabidopsis/10_45_45_all_samples/reads/Read1/spatial/merged_spatial/${sn}_filtered_R11.fastq
RV=Arabidopsis/10_45_45_all_samples/reads/Read2/tso_trimmed/${sn}_trimmed_filtered_R2.fastq

# References for mapping, annotation and ribo-filtering 
MAP=Genomes_transcriptomes_indexies/Arabidopsis/Arabidopsis_genome_mCitrine/StarIndex_Arabidopsis_Pseudomonas
ANN=Genomes_transcriptomes_indexies/Arabidopsis/Arabidopsis_genome_mCitrine/TAIR10_genes.gtf

# Barcodes settings 
ID=Genomes_transcriptomes_indexies/omni-v1_coordinates.txt

# Output folder and experiment name
OUTPUT=Arabidopsis/10_45_45_all_samples/data/${sn}
EXP=${sn}
TEMP=Arabidopsis/10_45_45_all_samples/temp/${sn}

# Running the pipeline

st_pipeline_run.py \
  --ids $ID \
  --ref-map $MAP \
  --expName $EXP \
  --ref-annotation $ANN \
  --mapping-threads 4 \
  --disable-multimap \
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
  --temp-folder $TEMP \
  --no-clean-up \
  --log-file $OUTPUT/${EXP}_log.txt \
  $FW $RV
