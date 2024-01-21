#!/bin/bash

############## Processing of RNA-Seq data
############## This file processes bulk RNA-seq samples starting from FASTQ to analysis-ready BAM output
############## It is sectioned into 5 parts: 0) Download RNA-Seq samples from NCBI (optional= if data is not downloaded) 1) Quality Check 2) Adapter trimming/Filtering Low Quality Read End Bases
############## 3) Mapping/Alignment 4) Mapping Quality Check
############## Input: RNA-Seq FASTQ file
############## Output: Analysis-ready BAM file

## Activate conda environment
module load anaconda3
source activate agshin_env

## Set required variables
export WD=/home/abalay/scratch/PD1I_datasets/

## Do not change these parameters unless they do not exist (Should be on server AGS_AB/Pipelines_and_Files/Packages or AGS_AB/Pipelines_and_Files/Files_for_scripts
export script_dir=/home/abalay/scratch/Pipelines_and_Files/Files_for_scripts


############## Part 0: Download RNA-Seq FASTQ samples from NCBI
Rscript $script_dir/RNASeq_processing_scripts/download_from_geo.r --parent-dir $WD --file-dir $WD/SRR_Acc_List.txt
$WD/downloadGEOscripts/download.jobsub.bat


############### Part 1: Quality Check
############### Input: FASTQ file
############### Output: FASTQC report files and MultiQC report html

## Quality Control (QC) of gunzip FASTQ file

### Run QC on FASTQ files
Rscript $script_dir/RNASeq_processing_scripts/RNASeq_fastqc.r --parent-dir $WD/raw  --output-dir $WD/processed
$WD/processed/scripts/fastqc/fastqc.jobsub.bat


### Summarize RNA-SEQ QC outcomes
### Create directory for multiqc reports and run multiqc for RNA-Seq QC
mkdir $WD/processed/multiqc
multiqc $WD/processed/sample_out/*/*_fastqc.zip --filename multiqc_qc_report --outdir $WD/processed/multiqc




############### Part 2: Adapter trimming & low quality base filtering (QUAL < 15, read maxN = 20 ) & short read removal (< 35 bp)
############### Input: FASTQ file
############### Output: Cutadapt processed FASTQ files & reports + MultiQC report html
Rscript $script_dir/RNASeq_processing_scripts/RNASeq_cutadapt.r --parent-dir $WD/raw/GSE213145  --output-dir $WD/processed
Rscript $script_dir/RNASeq_processing_scripts/RNASeq_cutadapt.r --parent-dir $WD/raw/GSE78220  --output-dir $WD/processed
Rscript $script_dir/RNASeq_processing_scripts/RNASeq_cutadapt.r --parent-dir $WD/raw/GSE91061  --output-dir $WD/processed

$WD/processed/scripts/cutadapt/cut.jobsub.bat

### Summarize RNA-SEQ Cutadapt outcomes
### Create directory for multiqc reports and run multiqc for RNA-Seq Cutadapt
multiqc $WD/processed/sample_out/*/*_fastqc.zip $WD/processed/sample_out/*/*_fastqc.zip $WD/processed/scripts/cutadapt/* --filename multiqc_cutadapt_report --outdir $WD/processed/multiqc




############### Part 3: Quality check (Part 2) to check adapter trimming and post-Cutadapt FASTQ read quality
############### Input: Cutadapt processed FASTQ file
############### Output: POST-Cutadapt FASTQC report files and MultiQC report html
Rscript $script_dir/RNASeq_processing_scripts/RNASeq_fastqc_postcutadapt.r --parent-dir $WD/processed
$WD/processed/scripts/refastqc/refastqc.jobsub.bat

### Summarize RNA-SEQ POST-Cutadapt FASTQC outcomes
### Create directory for multiqc reports and run multiqc for RNA-Seq post-Cutadapt QC
multiqc $WD/processed/sample_out/*/*_cutadapt_fastqc.zip --filename multiqc_cutadaptqc_report --outdir $WD/processed/multiqc/





############### Part 4: Genome alignment + statistics
############### Input: Cutadapt processed FASTQ file
############### Output: GSNAP generated BAM file & alignment statistics + MultiQC report html
Rscript $script_dir/RNASeq_processing_scripts/RNASeq_gsnap.r --parent-dir $WD/processed --splice-file $script_dir/human_splicesite_file.iit
$WD/processed/scripts/gsnap/gsnap.jobsub.bat

### Summarize RNA-SEQ GSNAP outcomes
### Create directory for multiqc reports and run multiqc for RNA-Seq GSNAP alignment
multiqc $WD/processed/sample_out/*/*.bam $WD/processed/sample_out/*/*.bam.bai $WD/processed/sample_out/*/*stats.txt --filename multiqc_align_report --outdir $WD/processed/multiqc/




############### Part 5: Alignment QC
############### Input: GSNAP processed BAM file
############### Output: Several reports with alignment statistics, read distribution, library strandedness & rRNA %
Rscript $script_dir/RNASeq_processing_scripts/RNASeq_rseqc.r --parent-dir $WD/processed --genome-anno-file /home/abalay/scratch/Pipelines_and_Files/Files_for_scripts/gencode.v37.annotation.bed --rrna-anno-file /home/abalay/scratch/Pipelines_and_Files/Files_for_scripts/hg38_rRNA.bed
$WD/processed/scripts/rseqc/rseqc.jobsub.bat

### Summarize RNA-SEQ alignment QC outcomes
### Create directory for multiqc reports and run multiqc for RNA-Seq alignment QC
multiqc $WD/processed/sample_out/*/*flagstats.txt $WD/processed/sample_out/*/*infstrand.txt $WD/processed/sample_out/*/*readistr.txt $WD/processed/sample_out/*/*stats.txt --filename multiqc_mappingqc_report --outdir $WD/processed/multiqc/


### Identify % of rRNA reads in each sample
chmod +x $script_dir/RNASeq_processing_scripts/rRNAcontamination.r
Rscript $script_dir/RNASeq_processing_scripts/rRNAcontamination.


