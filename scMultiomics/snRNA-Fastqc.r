#######################################
##################################
### FASTQC: quality check of snRNA-seq data for snMultiome samples (snATAC+snRNA-seq)
module load FastQC/0.11.7-Java-1.8.0_92
module load miniconda/4.12.0
# source activate /gpfs/gibbs/pi/kaminski/ab3478/conda_envs/agshin_R4env
source activate agshin_R4env

R

## Number of cores to multithread
num.cores <- 8
args <- commandArgs(trailingOnly=TRUE)

args[1] <- "/vast/palmer/scratch/kaminski/ab3478/snMultiome/snMultiome_RNA"
args[2] <- "/vast/palmer/scratch/kaminski/ab3478/snMultiome_analysis/snMultiome_RNA_analysis"

### Return error and cancel immediately if input directory doesn't exist
input.parent.dir <- args[1]
if(file.exists(input.parent.dir)==F){
    cat("Error: Input parent directory does not exist... Try again.\n")
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

fastqc.script.dir <- file.path(script.dir, "fastqc")
  if(file.exists(fastqc.script.dir)==F){
      cat("Generating fastqc script directory...\n")
      dir.create(fastqc.script.dir)
}

## wipe and remake the jobsub file
jobsub.filepath<-file.path(fastqc.script.dir,"fastqc.jobsub.bat")
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
        fastq.paths <-list.files(sample.paths[j], recursive=TRUE, full.names=TRUE)
#Return warning and skip sample if there are non-fastq files included in the folder
#### Distinguish Index1, Read1, Read2 and Read3
        index1.fastq.paths <-fastq.paths[grep("_I1_001.fastq.gz", fastq.paths)]
        read1.fastq.paths <-fastq.paths[grep("_R1_001.fastq.gz", fastq.paths)]
        read2.fastq.paths <-fastq.paths[grep("_R2_001.fastq.gz", fastq.paths)]
################################################# generate zUMI script for each sample
        index1.fastqs <- paste(index1.fastq.paths, sep = "", collapse = " ")
        read1.fastqs <- paste(read1.fastq.paths, sep = "", collapse = " ")
        read2.fastqs <- paste(read2.fastq.paths, sep = "", collapse = " ")
        sample.output.dir <- file.path(output.dir, sample.names[j])
        if(file.exists(sample.output.dir)==F){
            dir.create(sample.output.dir)
        }
        script.filepath <- file.path(fastqc.script.dir, paste(sample.names[j], "_fastqc.sh", sep = ""))
        sample.number <- paste0("S", j)
#### Make directory for merged I1, R1 and R2 files
        merged.dir <- file.path(sample.output.dir, "merged")
        if(file.exists(merged.dir)==F){
            dir.create(merged.dir)
        }
################## Handle I1 stuff
#### Make "I1 merged" directory
        # merged.I1.dir <- file.path(sample.output.dir, "merged_i1")
        # if(file.exists(merged.I1.dir)==F){
        #     dir.create(merged.I1.dir)
        # }
#### Make the a merged I1 file of I1.fastq.gz files from this sample in this batch
        if(length(index1.fastq.paths) > 1){
            cat("Hold on... merging Index1 files for ", sample.names[j], "...\n", sep = "")
        } else {
            cat("Hold on... transferring large Index1 file for ", sample.names[j], "...\n", sep = "")
        }
        merged.I1.file <- file.path(merged.dir, paste(sample.names[j], "merged", sample.number, "L003_I1_001.fastq.gz", sep = "_"))
######## Check if sample has replicate data already run w/in loop: if TRUE: append the merged fastq; if FALSE: overwrite
        if((sample.names[j]%in%completed.sample.names)==TRUE){
            cat("Appending the reads from sequencing replicates.\n")
            system(paste("cat ",index1.fastqs, " >> ", merged.I1.file, sep = ""))
        } else {
            system(paste("cat ",index1.fastqs, " > ", merged.I1.file, sep = ""))
        }
################## Handle R1 stuff
#### Make "R1 merged" directory
        # merged.R1.dir <- file.path(sample.output.dir, "merged_r1")
        # if(file.exists(merged.R1.dir)==F){
        #     dir.create(merged.R1.dir)
        # }
#### Make the a merged R1 file of R1.fastq.gz files from this sample in this batch
        if(length(read1.fastq.paths) > 1){
            cat("Hold on... merging Read1 files for ", sample.names[j], "...\n", sep = "")
        } else {
            cat("Hold on... transferring large Read1 file for ", sample.names[j], "...\n", sep = "")
        }
        merged.R1.file <- file.path(merged.dir, paste(sample.names[j], "merged", sample.number, "L003_R1_001.fastq.gz", sep = "_"))
######## Check if sample has replicate data already run w/in loop: if TRUE: append the merged fastq; if FALSE: overwrite
        if((sample.names[j]%in%completed.sample.names)==TRUE){
            cat("Appending the reads from sequencing replicates.\n")
            system(paste("cat ",read1.fastqs, " >> ", merged.R1.file, sep = ""))
        } else {
            system(paste("cat ",read1.fastqs, " > ", merged.R1.file, sep = ""))
        }
################# Handle R2 stuff
# #### Make "R2 merged" directory
#         merged.R2.dir <- file.path(sample.output.dir, "merged_r2")
#         if(file.exists(merged.R2.dir)==F){
#             dir.create(merged.R2.dir)
#         }
#### Make the a merged R2 file of R2.fastq.gz files from this sample in this batch
        if(length(read2.fastq.paths) > 1){
            cat("Hold on... merging Read2 files for ", sample.names[j], "...\n", sep = "")
        } else {
            cat("Hold on... transferring large Read2 file for ", sample.names[j], "...\n", sep = "")
        }
        merged.R2.file <- file.path(merged.dir, paste(sample.names[j], "merged", sample.number, "L003_R2_001.fastq.gz", sep = "_"))
######## Check if sample has replicate data already run w/in loop: if TRUE: append the merged fastq; if FALSE: overwrite
        if((sample.names[j]%in%completed.sample.names)==TRUE){
            cat("Appending the reads from sequencing replicates.\n")
            system(paste("cat ",read2.fastqs, " >> ", merged.R2.file, sep = ""))
        } else {
            system(paste("cat ",read2.fastqs, " > ", merged.R2.file, sep = ""))
        }
################################### Print the script
        cmd.out <- NULL
        cmd.out <- paste("#!/bin/bash\n\n")
        cmd.out <- paste(cmd.out,"#SBATCH --time=3:00:00 -p day,pi_kaminski --ntasks=1 --cpus-per-task=",num.cores,
            " --mem=51200M  --job-name=",paste(sample.names[j], "fastqc", sep = "."),
            " -o ",script.filepath,".o%J -e ",script.filepath,".e%J\n\n",sep="")
## Input
        cmd.out <- paste0(cmd.out, "fastqc ", read1.fastqs, " -o ", sample.output.dir, "\n")
        cmd.out <- paste0(cmd.out, "fastqc ", read2.fastqs, " -o ", sample.output.dir, "\n")
        cat(cmd.out,file=script.filepath,append=F)

#### Make sure the sample wasn't already processed as a replicate before appending sbatch to jobsub.bat
        if((sample.names[j]%in%completed.sample.names)==FALSE){
            cat("sbatch ",script.filepath,"\n",sep="",file=jobsub.filepath, append=TRUE)
        }
    ## record what's so we can check for replicate sequencing data runs and append their reads together
    completed.sample.names <- append(completed.sample.names, sample.names[j])
    }
}
system(paste("chmod 700 ",file.path(fastqc.script.dir,"fastqc.jobsub.bat"),sep=""))
#################
#completed.sample.names
cat("Completed generating .sh files for ", length(unique(completed.sample.names)), " samples.\n", sep="")
