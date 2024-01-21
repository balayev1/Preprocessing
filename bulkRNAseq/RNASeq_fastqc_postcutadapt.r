#######################################RNASeq_fastqc_postcutadapt.r
## argument 1: path to parent directory for the project
## e.g. Rscript RNASeq_fastqc_postcutadapt.r --parent-dir /home/abalay/scratch/PD1I_datasets/processed

##################################
### FASTQC: quality check of RNA-seq data post cutadapt

## Number of cores to multithread
num.cores <- 8

require("optparse")

### set the arguments
option_list = list(
    make_option(c("-p", "--parent-dir"), type="character", default=NULL, 
              help="parent directory with raw fastq files", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

### Return error and cancel immediately if input directory doesn't exist
input.parent.dir <- opt$p
if(file.exists(input.parent.dir)==F){
    cat("Error: Input parent directory does not exist... Try again.\n")
    exit()
}

sample.dir <- file.path(input.parent.dir, "sample_out")
if(file.exists(sample.dir)==F){
    cat("Generating sample directory...\n")
    dir.create(sample.dir)
}
script.dir <- file.path(input.parent.dir, "scripts")
if(file.exists(script.dir)==F){
    cat("Generating script directory...\n")
    dir.create(script.dir)
}

refastqc.script.dir <- file.path(script.dir, "refastqc")
  if(file.exists(refastqc.script.dir)==F){
      cat("Generating post cutadapt fastqc script directory...\n")
      dir.create(refastqc.script.dir)
}

## wipe and remake the jobsub file
jobsub.filepath<-file.path(refastqc.script.dir,"refastqc.jobsub.bat")
if(file.exists(jobsub.filepath)){
    cat("Creating new jobsub.bat...\n")
    unlink(jobsub.filepath)
}
file.create(jobsub.filepath)

## Extract the paths for the sample output directories
sample.names <- list.files(sample.dir)
sample.paths <- file.path(sample.dir, sample.names)

## Loop over samples
for(j in 1:length(sample.paths)){
#### return paths of all fastq files, recursively in case of subdirectories
#        fastq.paths <-list.files(file.path(sample.paths[j],raw.fastq.dir), full.names = TRUE, recursive = TRUE)
    fastq.paths <-list.files(sample.paths[j], recursive=TRUE, full.names=TRUE, pattern = "*_cutadapt.fastq.gz")
#### post cutadapt fastqc script for single-end RNA-seq data
    if (length(fastq.paths) == 1){
        sample.output.dir <- file.path(sample.dir, sample.names[j])
        if(file.exists(sample.output.dir)==F){
            dir.create(sample.output.dir)
        }
        script.filepath <- file.path(refastqc.script.dir, paste(sample.names[j], "_refastqc.sh", sep = ""))
################################### Print the script
        cmd.out <- NULL
        cmd.out <- paste("#!/bin/bash\n\n")
        cmd.out <- paste(cmd.out,"#SBATCH --time=3:00:00 --cpus-per-task=", num.cores,
        " --mem=16000M  --job-name=",paste(sample.names[j], "refastqc", sep = "."),
        " -o ",script.filepath,".o%J -e ",script.filepath,".e%J\n\n",sep="")
## Input
        cmd.out <- paste0(cmd.out, "fastqc ", fastq.paths, " -o ", sample.output.dir, "\n")
        cat(cmd.out,file=script.filepath,append=F)
        cat("sbatch ",script.filepath,"\n",sep="",file=jobsub.filepath, append=TRUE)
    }
#### post cutadapt fastqc script for paired-end RNA-seq data
    if (length(fastq.paths) > 1){
#### Distinguish Read1 and Read2
        read1.fastq.path <- fastq.paths[grep("_1_cutadapt.fastq.gz", fastq.paths)]
        read2.fastq.path <- fastq.paths[grep("_2_cutadapt.fastq.gz", fastq.paths)]

        sample.output.dir <- file.path(sample.dir, sample.names[j])
        if(file.exists(sample.output.dir)==F){
            dir.create(sample.output.dir)
        }
        script.filepath <- file.path(refastqc.script.dir, paste(sample.names[j], "_refastqc.sh", sep = ""))
        ################################### Print the script
        cmd.out <- NULL
        cmd.out <- paste("#!/bin/bash\n\n")
        cmd.out <- paste(cmd.out,"#SBATCH --time=3:00:00 --cpus-per-task=", num.cores,
        " --mem=16000M  --job-name=",paste(sample.names[j], "refastqc", sep = "."),
        " -o ",script.filepath,".o%J -e ",script.filepath,".e%J\n\n",sep="")
## Input
        cmd.out <- paste0(cmd.out, "fastqc ", read1.fastq.path, " -o ", sample.output.dir, "\n")
        cmd.out <- paste0(cmd.out, "fastqc ", read2.fastq.path, " -o ", sample.output.dir, "\n")
        cat(cmd.out,file=script.filepath,append=F)
        cat("sbatch ",script.filepath,"\n",sep="",file=jobsub.filepath, append=TRUE)
    }
}


system(paste("chmod 700 ",file.path(refastqc.script.dir,"refastqc.jobsub.bat"),sep=""))
