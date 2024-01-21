########################################################################RNASeq_cutadapt.R
## argument 1: parent directory for the relevant raw fastq files
## argument 2: parent directory for samples' respective output folders
## e.g. Rscript RNASeq_cutadapt.r --parent-dir /home/abalay/scratch/PD1I_datasets/raw --output-dir /home/abalay/scratch/PD1I_datasets/processed 
#####################################################################################################################################################
########## This Rscript is designed to use as input: a top level directory of RNA-Seq data
########## which contains 1 or more sub-directories, each representing a unique sequencing batch
########## which each contain sample specific subdirectories
########## which each contain raw fastq files
##########
##########
##########
########## Need to have your own cutadapt in PATH
##########
########## We filter out paired reads whose R1 or R2 was trimmed below 35bp length.
#####################################################################################################################################################


## Number of cores to multithread
num.cores <- 8

require("optparse")

### set the arguments
option_list = list(
    make_option(c("-p", "--parent-dir"), type="character", default=NULL, 
              help="parent directory with raw fastq files", metavar="character"),
    make_option(c("-o", "--output-dir"), type="character", default=NULL, 
              help="output directory for Cutadapt processed fastq files", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

### Return error and cancel immediately if input directory doesn't exist
input.parent.dir <- opt$p
if(file.exists(input.parent.dir)==F){
    cat("Error: Input parent directory does not exist... Try again.\n")
    exit()
}

## Generate the output directory structure
output.parent.dir <- opt$o
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
if(file.exists(jobsub.filepath) == FALSE){
    cat("Creating new jobsub.bat...\n")
    file.create(jobsub.filepath)
}

## Extract the path for the data parent directory (batches)
# soup.batch.names <- list.files(file.path(input.parent.dir))
soup.batch.dir.paths <- file.path(input.parent.dir)

## make blank record for processed samples outside of the loop
completed.sample.names <- vector()

## Loop over batches
for(i in 1:length(soup.batch.dir.paths)){
    sample.names <- list.files(soup.batch.dir.paths[i])
    sample.paths <- file.path(soup.batch.dir.paths[i], sample.names)
## Loop over samples
    for(j in 1:length(sample.paths)){
        #        fastq.paths <-list.files(file.path(sample.paths[j],raw.fastq.dir), full.names = TRUE, recursive = TRUE)
        fastq.paths <-list.files(sample.paths[j], recursive=TRUE, full.names=TRUE)
        #Return warning and skip sample if there are non-fastq files included in the folder
        #### cutadapt script for single-end RNA-seq data
        if (length(fastq.paths) == 1){
            sample.output.dir <- file.path(output.dir, sample.names[j])
            if(file.exists(sample.output.dir)==F){
                dir.create(sample.output.dir)
            }
            script.filepath <- file.path(cut.script.dir, paste(sample.names[j], "_cutadapt.sh", sep = ""))
        ################################### Print the script
            cmd.out <- NULL
            cmd.out <- paste("#!/bin/bash\n\n")
            cmd.out <- paste(cmd.out,"#SBATCH --time=5:00:00 --cpus-per-task=", num.cores,
            " --mem=16000M  --job-name=",paste(sample.names[j], "cutadapt", sep = "."),
            " -o ",script.filepath,".o%J -e ",script.filepath,".e%J\n\n",sep="")
        ## Input
            cmd.out <- paste0(cmd.out, "cutadapt -a AGATCGGAAGAG -q 15 -m 35 --max-n 20 -j ", 
                num.cores, 
                " -o ", file.path(sample.output.dir, paste0(sample.names[j], "_cutadapt.fastq.gz ")), fastq.paths, "\n")
            cat(cmd.out,file=script.filepath,append=F)
        #### Make sure the sample wasn't already processed as a replicate before appending sbatch to jobsub.bat
            if((sample.names[j]%in%completed.sample.names)==FALSE){
                cat("sbatch ",script.filepath,"\n",sep="",file=jobsub.filepath, append=TRUE)
            }
        ## record what's so we can check for replicate sequencing data runs and append their reads together
            completed.sample.names <- append(completed.sample.names, sample.names[j])
        }
#### cutadapt script for paired-end RNA-seq data
        if (length(fastq.paths) > 1){
#### Distinguish Read1 and Read2
            read1.fastq.path <- fastq.paths[grep("_1.fastq", fastq.paths)]
            read2.fastq.path <- fastq.paths[grep("_2.fastq", fastq.paths)]

            sample.output.dir <- file.path(output.dir, sample.names[j])
            if(file.exists(sample.output.dir)==F){
                dir.create(sample.output.dir)
            }
            script.filepath <- file.path(cut.script.dir, paste(sample.names[j], "_cutadapt.sh", sep = ""))
            ################################### Print the script
            cmd.out <- NULL
            cmd.out <- paste("#!/bin/bash\n\n")
            cmd.out <- paste(cmd.out,"#SBATCH --time=5:00:00 --cpus-per-task=", num.cores,
            " --mem=16000M  --job-name=",paste(sample.names[j], "cutadapt", sep = "."),
            " -o ",script.filepath,".o%J -e ",script.filepath,".e%J\n\n",sep="")
## Input
            cmd.out <- paste0(cmd.out, "cutadapt -a AGATCGGAAGAG -q 15 -m 35 --max-n 20 --pair-filter=any -j ", num.cores,
                " -o ", file.path(sample.output.dir, paste0(sample.names[j], "_1_cutadapt.fastq.gz ")), 
                " -p ", file.path(sample.output.dir, paste0(sample.names[j], "_2_cutadapt.fastq.gz ")), " ", 
                read1.fastq.path, " ", read2.fastq.path, "\n")
            cat(cmd.out,file=script.filepath,append=F)
#### Make sure the sample wasn't already processed as a replicate before appending sbatch to jobsub.bat
            if((sample.names[j]%in%completed.sample.names)==FALSE){
                cat("sbatch ",script.filepath,"\n",sep="",file=jobsub.filepath, append=TRUE)
            }
## record what's so we can check for replicate sequencing data runs and append their reads together
            completed.sample.names <- append(completed.sample.names, sample.names[j])
        }
    }
}  
system(paste("chmod 700 ",file.path(cut.script.dir,"cut.jobsub.bat"),sep=""))
#################
#completed.sample.names
cat("Completed generating .sh files for ", length(unique(completed.sample.names)), " samples.\n", sep="")