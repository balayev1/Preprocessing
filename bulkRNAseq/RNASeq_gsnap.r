## argument 1: path to parent directory for the project
## eg: Rscript RNASeq_gsnap.r --parent-dir $WD/processed --splice-file /home/abalay/scratch/Pipelines_and_Files/Files_for_scripts/human_splicesite_file.iit
#####################################################################################################################################################
########## This Rscript is designed to use as input: the same directory labeled as "output.parent.dir" in the cutadapt Rscript
########## which contains a directory titles "sample_out" and "scripts"
########## The cutadapt job must be run beforehand to generate cutadapt output fastq files in the
##########
########## The purpose of the script is to generate a sample-specific .sh for processing each sample through GSNAP
##########
########## !!! Note: Should have GSNAP (GMAP aligner) installed
########## !!! Must generate splice site file and genome index (See GMAP aligner tutorial webpage)
########## GRCh38.p13 Human genome used
##########
###### Summary:
###### 1 Identify all sample folders
#loop{ 2 generate script for GSNAP}
#####################################################################################################################################################
# 8 threads
# 60GB memory

num.cores <- 8

require("optparse")

### set the arguments
option_list = list(
    make_option(c("-p", "--parent-dir"), type="character", default=NULL, 
              help="parent directory with raw fastq files", metavar="character"),
    make_option(c("-s", "--splice-file"), type="character", default=NULL, 
              help="directory to splice site file from GSNAP (must be in iit format see GMAP aligner instructions)", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

## set parent directory
parent.dir <- opt$p

## set sample and script output directories
sample.parent.dir <- file.path(parent.dir, "sample_out")
scripts.parent.dir <- file.path(parent.dir, "scripts")

gsnap.scripts.dir <- file.path(scripts.parent.dir, "gsnap")
if(file.exists(gsnap.scripts.dir)==F){
    dir.create(gsnap.scripts.dir)
}

soup.sample.names <- list.files(file.path(sample.parent.dir))
soup.sample.dir.paths <- file.path(sample.parent.dir, soup.sample.names)
jobsub.filepath <- file.path(gsnap.scripts.dir,"gsnap.jobsub.bat")

gsnap.path <- "gsnap"

splice.site.path <- opt$s
# anno.path <- "/home/abalay/scratch/Pipelines_and_Files/Files_for_scripts/human_annotation.iit"

gsnap.options <- paste("--gunzip -d GRCh38 --dir /home/abalay/data/conda/envs/agshin_env/share/",
    " --nthreads ", num.cores, " --novelsplicing 1", sep="") # enable novel splice site option 

## Wipe the old jobsub.bat file
cat("",sep="",file=jobsub.filepath, append=FALSE)
## Loop over samples to extract path to cutadapt processed files
for(i in 1:length(soup.sample.names)){
    cutadapt.fastq.files <- list.files(soup.sample.dir.paths[i], pattern = "*_cutadapt.fastq.gz")
    cutadapt.fastq.paths <- file.path(soup.sample.dir.paths[i], cutadapt.fastq.files)

#### gsnap script for single-end RNA-seq data
    if (length(cutadapt.fastq.paths) == 1){
        script.filepath <- file.path(gsnap.scripts.dir, paste(soup.sample.names[i], "_gsnap.sh", sep = ""))
        ##### print the script -g stands for genome annotation;
        cmd.out <- NULL
        cmd.out <- paste("#!/bin/bash\n")
        cmd.out <- paste(cmd.out,"#SBATCH --time=48:00:00 --cpus-per-task=", num.cores,
            " --mem=60000M --job-name=", paste0(soup.sample.names[i], "gsnap"),
            " -o ",script.filepath,".o%J -e ",script.filepath,".e%J\n",sep="")
        cmd.out <- paste(cmd.out, gsnap.path, " ", gsnap.options, " --use-splicing ", splice.site.path, 
            " --format sam ", cutadapt.fastq.paths, "| samtools view -hb| samtools sort -@ ", num.cores,
            " > ", file.path(soup.sample.dir.paths[i],paste(soup.sample.names[i], "_cutadapt_sorted.bam\n", sep="")), sep="")
        cmd.out <- paste(cmd.out, "samtools index ", 
            file.path(soup.sample.dir.paths[i], paste0(soup.sample.names[i], "_cutadapt_sorted.bam")), " > ", 
            file.path(soup.sample.dir.paths[i],paste0(soup.sample.names[i], "_cutadapt_sorted.bam.bai\n")), sep="")
        cmd.out <- paste(cmd.out, "samtools stats ", 
            file.path(soup.sample.dir.paths[i], paste0(soup.sample.names[i], "_cutadapt_sorted.bam")), " > ",
            file.path(soup.sample.dir.paths[i],paste0(soup.sample.names[i], "_cutadapt_sorted.stats.txt\n")), sep="")
        cat(cmd.out,file=script.filepath,append=F)
        cat("sbatch ",script.filepath,"\n",sep="",file=jobsub.filepath, append=TRUE)
    }
#### gsnap script for paired-end RNA-seq data
    if (length(cutadapt.fastq.paths) > 1){
        script.filepath <- file.path(gsnap.scripts.dir, paste(soup.sample.names[i], "_gsnap.sh", sep = ""))
        #### Distinguish Read1 and Read2
        read1.cut.fastq.path <- cutadapt.fastq.paths[grep("_1_cutadapt.fastq.gz", cutadapt.fastq.paths)]
        read2.cut.fastq.path <- cutadapt.fastq.paths[grep("_2_cutadapt.fastq.gz", cutadapt.fastq.paths)]

        ##### print the script -g stands for genome annotation;
        cmd.out <- NULL
        cmd.out <- paste("#!/bin/bash\n")
        cmd.out <- paste(cmd.out,"#SBATCH --time=48:00:00 --cpus-per-task=", num.cores,
            " --mem=60000M --job-name=", paste0(soup.sample.names[i], "gsnap"),
            " -o ",script.filepath,".o%J -e ",script.filepath,".e%J\n",sep="")
        cmd.out <- paste(cmd.out, gsnap.path, " ", gsnap.options, " --use-splicing ", splice.site.path, 
            " --format sam ", read1.cut.fastq.path, " ", read2.cut.fastq.path, " ",
            "| samtools view -hb| samtools sort -@ ", num.cores, " > ", 
            file.path(soup.sample.dir.paths[i],paste0(soup.sample.names[i], "_cutadapt_sorted.bam\n")), sep="")
        cmd.out <- paste(cmd.out, "samtools index ", 
            file.path(soup.sample.dir.paths[i], paste0(soup.sample.names[i], "_cutadapt_sorted.bam")), " > ", 
            file.path(soup.sample.dir.paths[i],paste0(soup.sample.names[i], "_cutadapt_sorted.bam.bai\n")), sep="")
        cmd.out <- paste(cmd.out, "samtools stats ", 
            file.path(soup.sample.dir.paths[i], paste0(soup.sample.names[i], "_cutadapt_sorted.bam")), " > ",
            file.path(soup.sample.dir.paths[i],paste0(soup.sample.names[i], "_cutadapt_sorted.stats.txt\n")), sep="")
        cat(cmd.out,file=script.filepath,append=F)
        cat("sbatch ",script.filepath,"\n",sep="",file=jobsub.filepath, append=TRUE)
    }
}
system(paste("chmod 700 ",jobsub.filepath,sep=""))