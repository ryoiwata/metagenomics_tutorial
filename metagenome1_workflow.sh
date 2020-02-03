#!/bin/bash

# USAGE: sh metagenome1_workflow.sh <fastq file> 
# This script accepts the location of a fastq file and performs the following steps on it:
    ## starting with trimming
    ## convert the  

# #start a conda environment
source ~/anaconda3/etc/profile.d/conda.sh
conda env list

# load the modules
module load bbtools/37.02
module load megahit/1.1.1
module load samtools/1.4.1
module load quast/3.1

# Get input file and locations  
input_original_file_name=$1
input_file_path=$(realpath $input_original_file_name)
echo "Input file's path: $input_file_path"

# Grab base of filename for future naming
input_basename=$(basename -- $input_file_path)
echo "Starting analysis of: $input_basename"

#Creaing an array of the basename split by periods
IFS='.' read -r -a array <<< "$input_basename"

#Remove extension
let "basename_index = ${#array[@]} -2"
basename_no_ext=$(echo $input_basename | cut -d "." -f 1-$basename_index)

# Changing to our working directory
project_path=$(pwd)
mkdir $project_path/metagenome1
cd $project_path/metagenome1
mkdir $project_path/metagenome1/data
mkdir $project_path/metagenome1/data/stdout

echo "Step 1: Filtering contaminants"
# Filter out contaminant reads placing them in their own file
bbduk.sh in=${input_file_path} \
out=${project_path}/metagenome1/data/${basename_no_ext}.unmatched.fq.gz \
outm=${project_path}/metagenome1/data/${basename_no_ext}.matched.fq.gz \
k=31 \
hdist=1 \
ftm=5 \
ref=/software/apps/bbtools/gcc/64/37.02/resources/sequencing_artifacts.fa.gz \
stats=${project_path}/metagenome1/data/contam_stats.txt \
2>&1 | tee $project_path/metagenome1/data/stdout/contaminant_filter.txt

echo "Step 2: Trimming Adapters"
# Trim the adapters using the reference file adaptors.fa (provided by bbduk)
bbduk.sh in=${project_path}/metagenome1/data/${basename_no_ext}.unmatched.fq.gz \
out=${project_path}/metagenome1/data/${basename_no_ext}.trimmed.fastq.gz \
ktrim=r \
k=23 \
mink=11 \
hdist=1 \
tbo=t \
qtrim=r \
trimq=20 \
ref=/software/apps/bbtools/gcc/64/37.02/resources/adapters.fa \
2>&1 | tee $project_path/metagenome1/data/stdout/adapter_trim.txt

echo "Step 3: Merging Reads"
# Merge the reads together
bbmerge.sh in=${project_path}/metagenome1/data/${basename_no_ext}.trimmed.fastq.gz \
out=${project_path}/metagenome1/data/${basename_no_ext}.merged.fastq.gz \
outu=${project_path}/metagenome1/data/${basename_no_ext}.unmerged.fastq.gz \
ihist=${project_path}/metagenome1/data/insert_size.txt \
usejni=t \
2>&1 | tee $project_path/metagenome1/data/stdout/merge_reads.txt

echo "Step 4: Metagenome Assembly"
#metagenome assembly
megahit  -m 1200000000000 \
-t 38 \
--read ${project_path}/metagenome1/data/${basename_no_ext}.merged.fastq.gz \
--12 ${project_path}/metagenome1/data/${basename_no_ext}.unmerged.fastq.gz \
--k-list 21,41,61,81,99 \
--no-mercy \
--min-count 2 \
--out-dir ${project_path}/metagenome1/data/megahit/ \
2>&1 | tee $project_path/metagenome1/data/stdout/metagenome_assembly.txt


#fastg file for vizualization
megahit_toolkit contig2fastg 99 \
${project_path}/metagenome1/data/megahit/intermediate_contigs/k99.contigs.fa > ${project_path}/metagenome1/data/k99.fastg

echo "Step 5: Mapping Reads to a Reference"
#map reads back to known reference
bbmap.sh \
  ref=/project/microbiome_workshop/metagenome/viral/data/Ref_database.fna \
  in=${project_path}/metagenome1/data/${basename_no_ext}.trimmed.fastq.gz \
  out=${project_path}/metagenome1/data/${basename_no_ext}.map_to_genomes.sam \
  nodisk \
  usejni=t \
  pigz=t \
  threads=40 \
  bamscript=bamscript1.sh \
  2>&1 | tee $project_path/metagenome1/data/stdout/map_read_reference.txt

#Convert SAM file to BAM file
${project_path}/metagenome1/bamscript1.sh

echo "Step 6: Mapping Reads to Contigs"
#map reads back to contigs
bbmap.sh \
  ref=data/megahit/final.contigs.fa \
  in=${project_path}/metagenome1/data/${basename_no_ext}.trimmed.fastq.gz \
  out=${project_path}/metagenome1/data/${basename_no_ext}.map_to_contigs.sam \
  nodisk \
  usejni=t \
  pigz=t \
  threads=40 \
  bamscript=bamscript2.sh \
  2>&1 | tee $project_path/metagenome1/data/stdout/map_read_contig.txt


#Convert SAM file to BAM file
${project_path}/metagenome1/bamscript2.sh

echo "job complete!"