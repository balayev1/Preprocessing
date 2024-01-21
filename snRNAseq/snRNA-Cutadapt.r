########################################################################10X_cutadapt.R
## argument 1: parent directory for the relevant raw fastq files
## argument 2: parent directory for samples' respective output folders
## e.g. Rscript 10X_cutadapt.R /gpfs/loomis/project/fas/kaminski/public/Backup/Taylor/10X/CellSoups /home/fas/kaminski/tsa27/scratch60/10x_zUMI
#####################################################################################################################################################
########## This Rscript is designed to use as input: a top level directory of 10X visium data
########## which contains 1 or more sub-directories, each representing a unique sequencing batch
########## which each contain sample specific subdirectories
########## which each contain raw fastq taken from bcl2fastq (or cellranger mkfastq); often, but not always within a subfolder titled "raw_fastqs"
##########
##########
##########
########## Need to have your own cutadapt in PATH
########## I'm using cutatdapt 1.17 (installation code below)
##########
########## We filter out paired reads whose R2 was trimmed below 30bp length.
########## Output is prepared to be passed to the zUMIs pipeline.
##########
###### Summary:
###### 1 generate output file structure for pipeline
#loop{ 2 merge multiple R1 and R2 files into individual files
#loop# 3 merge R1 and R2 files from replicate library sequencing runs
#loop} 4 generate scripts (& batch executable) to run each samples' merged R1 & R2 files through cutadapt
#####################################################################################################################################################
## Install cutadapt 1.17 in local PATH (can't get a functional one in a public directory)
#/gpfs/loomis/project/fas/kaminski/public/softwares/Python-3.6.6/bin/pip3 install --user  cutadapt
## using cutadapt version 1.17
## activate conda environment
module load miniconda/4.12.0
# source activate /gpfs/gibbs/pi/kaminski/ab3478/conda_envs/agshin_R4env
source activate agshin_R4env

R

## Number of cores to multithread
num.cores <- 8
args <- commandArgs(trailingOnly=TRUE)

args[1] <- "/vast/palmer/scratch/kaminski/ab3478/Covid19_Haberman_2020/rawData"
args[2] <- "/vast/palmer/scratch/kaminski/ab3478/Covid19_Haberman_2020/processedData"

### Return error and cancel immediately if input directory doesn't exist
input.parent.dir <- args[1]
if(file.exists(input.parent.dir)==F){
    cat("Error: Input parent directory does not exist... Try again.\n")
    exit()
}

## Generate the output directory structure
output.parent.dir <- args[2]
if(file.exists(output.parent.dir)==F){
    cat("Generating top level directory for output...\n")
    dir.create(output.parent.dir)
}
output.dir <- file.path(output.parent.dir, "sample_out")
if(file.exists(output.dir)==F){
    cat("Generating parent directory for sample output...\n")
    dir.create(output.dir)
}
script.dir <- file.path(output.parent.dir, "scripts")
if(file.exists(script.dir)==F){
    cat("Generating script directory...\n")
    dir.create(script.dir)
}

cut.script.dir <- file.path(script.dir, "cutadapt")
  if(file.exists(cut.script.dir)==F){
      cat("Generating cutadapt script directory...\n")
      dir.create(cut.script.dir)
}

## wipe and remake the jobsub file
jobsub.filepath<-file.path(cut.script.dir,"cut.jobsub.bat")
if(file.exists(jobsub.filepath)){
    cat("Creating new jobsub.bat...\n")
    file.remove(jobsub.filepath)
}
file.create(jobsub.filepath)

##Extract the paths for the 10x CellSoup data parent directories (batches)
soup.batch.names <- list.files(file.path(input.parent.dir))
soup.batch.dir.paths <- file.path(input.parent.dir, soup.batch.names)
## make blank record for processed samples outside of the loop
completed.sample.names <- vector()
## Loop over batches
for(i in 1:length(soup.batch.dir.paths)){
    sample.names <- list.files(soup.batch.dir.paths[i])
    sample.paths <- file.path(soup.batch.dir.paths[i], sample.names)
## Loop over samples
    for(j in 1:length(sample.paths)){
#### return paths of all fastq files, recursively in case of subdirectories
#        fastq.paths <-list.files(file.path(sample.paths[j],raw.fastq.dir), full.names = TRUE, recursive = TRUE)
        fastq.paths <-list.files(sample.paths[j], recursive=TRUE, full.names=TRUE)
#Return warning and skip sample if there are non-fastq files included in the folder
#### Distinguish Read1 and Read2
        read1.fastq.paths <-fastq.paths[grep("_R1_", fastq.paths)]
        read2.fastq.paths <-fastq.paths[grep("_R2_", fastq.paths)]
################################################# generate zUMI script for each sample
        read1.fastqs <- paste(read1.fastq.paths, sep = "", collapse = " ")
        read2.fastqs <- paste(read2.fastq.paths, sep = "", collapse = " ")
        sample.output.dir <- file.path(output.dir, sample.names[j])
        if(file.exists(sample.output.dir)==F){
            dir.create(sample.output.dir)
        }
        script.filepath <- file.path(cut.script.dir, paste(sample.names[j], "_cutadapt.sh", sep = ""))
################## Handle R1 stuff
#### Make "R1 merged" directory
        merged.R1.dir <- file.path(sample.output.dir, "merged_r1")
        if(file.exists(merged.R1.dir)==F){
            dir.create(merged.R1.dir)
        }
#### Make the a merged R1 file of R1.fastq.gz files from this sample in this batch
        merged.R1.file <- file.path(merged.R1.dir, paste(sample.names[j], "_R1_merged.fastq.gz", sep = ""))
################# Handle R2 stuff
#### Make "R2 merged" directory
        merged.R2.dir <- file.path(sample.output.dir, "merged_r2")
        if(file.exists(merged.R2.dir)==F){
            dir.create(merged.R2.dir)
        }
#### Make the a merged R2 file of R2.fastq.gz files from this sample in this batch
        merged.R2.file <- file.path(merged.R2.dir, paste(sample.names[j], "_R2_merged.fastq.gz", sep = ""))
################################### Make the output directory for merged and trimmed fastqs
        trimmed.merged.fastq.dir <- file.path(sample.output.dir,"trimmed_merged_fastq")
        if(file.exists(trimmed.merged.fastq.dir)==F){
            dir.create(trimmed.merged.fastq.dir)
        }
################################### Print the script
        cmd.out <- NULL
        cmd.out <- paste("#!/bin/bash\n\n")
        cmd.out <- paste(cmd.out,"#SBATCH --time=3:00:00 -p day,pi_kaminski,bigmem --ntasks=1 --cpus-per-task=",num.cores,
        " --mem=51200M  --job-name=",paste(sample.names[j], "cut", sep = "."),
        " -o ",script.filepath,".o%J -e ",script.filepath,".e%J\n\n",sep="")
## Trim TSO reverse complement (SMART II A oligo); interleave the output and pipe it to another cutadapt run
        cmd.out<- paste0(cmd.out, "cutadapt -G AAGCAGTGGTATCAACGCAGAGTACTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT --cores=",num.cores," --interleaved ",
## Input
        merged.R1.file," ",
        merged.R2.file, " | \n\n",
####### Run cutadapt again; trim poly A (up to 2 errors); trim poly G (up to 3 errors); discard reads with >20 N; discard reads w/ R2 length < 30bp
        "cutadapt -G \"T{20}\" -G \"G{30}\" --max-n 20 --pair-filter=any --minimum-length=30 --cores=",num.cores," --interleaved",
## Output
        " -o ",file.path(trimmed.merged.fastq.dir, "trimmed.R1.fastq.gz"),
        " -p ",file.path(trimmed.merged.fastq.dir, "trimmed.R2.fastq.gz -"))
    cat(cmd.out,file=script.filepath,append=F)
#### Make sure the sample wasn't already processed as a replicate before appending sbatch to jobsub.bat
    if((sample.names[j]%in%completed.sample.names)==FALSE){
        cat("sbatch ",script.filepath,"\n",sep="",file=jobsub.filepath, append=TRUE)
    }
## record what's so we can check for replicate sequencing data runs and append their reads together
   completed.sample.names <- append(completed.sample.names, sample.names[j])
    }
}  
system(paste("chmod 700 ",file.path(cut.script.dir,"cut.jobsub.bat"),sep=""))
#################
#completed.sample.names
cat("Completed generating .sh files for ", length(unique(completed.sample.names)), " samples.\n", sep="")

