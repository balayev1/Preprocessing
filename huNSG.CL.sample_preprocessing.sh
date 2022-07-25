#!/bin/bash

############## Preprocessing of huNSG mice-derived cell lines and related infected cell lines with EBV and/or KSHV


############## This file processes bulk RNA-seq samples starting from FASTQ to analysis-ready BAM output
############## It is sectioned into 5 parts: 0) Download RNA-Seq samples from NCBI (optional= if data is not downloaded) 1) Quality Check 2) Adapter trimming/Filtering Low Quality Read End Bases
############## 3) Mapping/Alignment 4) Mapping Quality Check
############## Input: RNA-Seq FASTQ file
############## Output: Analysis-ready BAM file

## Load the environment 
module load hpc
srun --mem-per-cpu=4GB --time=6:00:00 --pty --cpus-per-task=4 bash 
module load anaconda3
source activate agshin_condaenv


## Set required variables
export WD=/home/cluster/abalay/scratch/Cell_Lines

## Do not change these parameters unless they do not exist (Should be on server AGS_AB/Pipelines_and_Files/Packages or AGS_AB/Pipelines_and_Files/Files_for_scripts
export script_dir=/home/cluster/abalay/scratch/Pipelines_and_Files
export gsnap=$script_dir/Packages/gmap-2020-06-01/src/gsnap
export ref_gen=$script_dir/Files_for_scripts/GRCh38.p13.genome.fa 
export rseqc=$script_dir/Packages/RSeQC-3.0.1/scripts
export gsnap=$script_dir/Packages/gmap-2020-06-01/src/gsnap

####### Make folder to keep scripts
mkdir $WD/Scripts


############## Part 0: Download RNA-Seq FASTQ samples from NCBI
python3 $WD/Scripts/file_download_generation.sh --file-path $WD/file_download_script.sh
chmod +x $WD/final_file_download.sh 
$WD/final_file_download.sh

# Check number of folder per each sample
ls -d */ | wc -l
#92-1=91




############### Part 1: Quality Check
############### Input: FASTQ file
############### Output: FASTQC report files and MultiQC report html

### Run QC on FASTQ files
python3 $WD/Scripts/fastqc.py --work-dir $WD --script-dir $script_dir --output-dir $WD
chmod +x $WD/schedule_jobs_fastqc.sh
$WD/schedule_jobs_fastqc.sh

# Check number of FASTQ reports
find $WD -name "*_fastqc.zip" |wc -l
#170

### Create directory for multiqc reports and run multiqc
mkdir $WD/multiqc_reports
multiqc -o  $WD/multiqc_reports/cell_lines_fastqc_step1 $(find $WD -name "*_fastqc.zip")







################# Part 2: Adapter trimming/Filtering Low Quality Read End Bases
################# Input: FASTQ file
################# Output: processed FASTQ and report file


## Adapter trimming & Filtering low-quality 3'-end read bases
### Make subdirectory for operation
python3 $WD/Scripts/cutadapt.py --work-dir $WD --output-dir $WD
chmod +x $WD/schedule_jobs_cutadapt.sh
$WD/schedule_jobs_cutadapt.sh

# Check number of cutadapt processed fastq files
find $WD -name "*_cutadapt_fastq.gz" | wc -l
# 170



### Check adapter trimming quality and overrepresented sequences
#### Make post-trimming FASTQC report
python3 $WD/Scripts/fastqc_post_cutadapt.py --work-dir $WD --script-dir $script_dir --output-dir $WD
chmod +x $WD/schedule_jobs_fastqc_post_cutadapt.sh
$WD/schedule_jobs_fastqc_post_cutadapt.sh



### Run MultiQC to observe report of all FASTQ files on one html page
multiqc -o  $WD/multiqc_reports/cell_lines_fastqc_step2 $(find $WD -name "*cutadapt_fastq_fastqc.zip") $(find $WD -name "*cutadapt.out") $(find $WD -name "*cutadapt.fastq.gz")









########### Part 3: Mapping/Alignment
########### Input: Processed FASTQ file
########### Output: BAM, BAM index
python3 $WD/Scripts/gsnap.py --work-dir $WD --script-dir $script_dir --output-dir $WD
chmod +x $WD/schedule_jobs_gsnap.sh
$WD/schedule_jobs_gsnap.sh

# Check number of BAM files
find $WD -name "*.bam" | wc -l
#91

# Check number of index BAM.BAI iles
find $WD -name "*.bam.bai" | wc -l
#91



################ Part 4: Mapping Quality Check
################ Input: BAM file
################ Output: Statistics reports for Mapping Quality Check
python3 $WD/Scripts/mapping_qc.py --work-dir $WD --script-dir $script_dir --output-dir $WD
chmod +x $WD/schedule_jobs_mapping_qc.sh
$WD/schedule_jobs_mapping_qc.sh


### Run MultiQC to observe report of all FASTQ files on one html page
multiqc -o  $WD/multiqc_reports/cell_lines_fastqc_step3 $(find $WD -name "*bam") $(find $WD -name "*txt")


### Check rRNA contamination
Rscript $WD/Scripts/rRNAcontamination.R


