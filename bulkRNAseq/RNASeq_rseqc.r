## argument 1: path to parent directory for the project
## eg: Rscript RNASeq_rseqc.r --parent-dir $WD/processed --genome-anno-file /home/abalay/scratch/Pipelines_and_Files/Files_for_scripts/gencode.v37.annotation.bed --rrna-anno-file /home/abalay/scratch/Pipelines_and_Files/Files_for_scripts/hg38_rRNA.bed
#####################################################################################################################################################
########## This Rscript is designed to use as input: the argument directory in the gsnap Rscript
########## which contains a directory titles "sample_out" and "scripts"
########## The gsnap job must be run beforehand to generate aligned output bam files in the sample directory "sample_out"
##########
########## The purpose of the script is to generate a sample-specific .sh for processing each sample through RSEQC
########## to check the alignment quality of the reads, infer the library strandedness (optional), check the read distribution
########## across genome and verify rRNA %

########## !!! Note: PACKAGES RSEQC and SAMTOOLS ARE NEEDED TO RUN THE SCRIPT

###### Summary:
###### 1 Identify all sample folders
#loop{ 2 generate script for RSEQC}
#####################################################################################################################################################
# 8 threads
# 60GB memory


num.cores <- 8

require("optparse")

### set the arguments
option_list = list(
    make_option(c("-p", "--parent-dir"), type="character", default=NULL, 
              help="parent directory with raw fastq files", metavar="character"),
    make_option(c("-g", "--genome-anno-file"), type="character", default=NULL, 
              help="parent directory with raw fastq files", metavar="character"),
    make_option(c("-r", "--rrna-anno-file"), type="character", default=NULL, 
              help="parent directory with raw fastq files", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

## set parent directory
parent.dir <- opt$p

sample.parent.dir <- file.path(parent.dir, "sample_out")
scripts.parent.dir <- file.path(parent.dir, "scripts")

rseqc.scripts.dir <- file.path(scripts.parent.dir, "rseqc")
if(file.exists(rseqc.scripts.dir)==F){
    dir.create(rseqc.scripts.dir)
}

soup.sample.names <- list.files(file.path(sample.parent.dir))
soup.sample.dir.paths <- file.path(sample.parent.dir, soup.sample.names)
jobsub.filepath <- file.path(rseqc.scripts.dir,"rseqc.jobsub.bat")

# rseqc.path <- "rseqc"
samtools.path <- "samtools"

genome.bedfile.path <- opt$g
hgrna.bedfile.path <- opt$r

## Wipe the old jobsub.bat file
cat("",sep="",file=jobsub.filepath, append=FALSE)
## Loop over samples to extract path to cutadapt processed files
for(i in 1:length(soup.sample.names)){
    bam.file <- list.files(soup.sample.dir.paths[i], pattern = "*_cutadapt_sorted.bam")[1]
    bam.file.path <- file.path(soup.sample.dir.paths[i], bam.file)[1]

#### rseqc script for RNA-seq data
    script.filepath <- file.path(rseqc.scripts.dir, paste(soup.sample.names[i], "_rseqc.sh", sep = ""))
    ##### print the script;
    cmd.out <- NULL
    cmd.out <- paste("#!/bin/bash\n")
    cmd.out <- paste(cmd.out,"#SBATCH --time=12:00:00 --cpus-per-task=", num.cores,
        " --mem=60000M --job-name=", paste0(soup.sample.names[i], "gsnap"),
        " -o ",script.filepath,".o%J -e ",script.filepath,".e%J\n",sep="")
    ### print flagstat script
    cmd.out <- paste(cmd.out, samtools.path, " flagstat --threads ", num.cores, " ", bam.file.path , " > ",
        file.path(soup.sample.dir.paths[i], paste0(strsplit(bam.file, ".bam")[[1]], ".flagstats.txt\n")), sep="")
    
    ### print infer strandedness script
    cmd.out <- paste(cmd.out, "infer_experiment.py -r ", genome.bedfile.path, 
        " -s 2000000  -q 10 -i ", bam.file.path, " > ", 
        file.path(soup.sample.dir.paths[i], paste0(strsplit(bam.file, ".bam")[[1]], ".infstrand.txt\n")), sep="")

    ### print read distribution script
    cmd.out <- paste(cmd.out, "read_distribution.py -r ", genome.bedfile.path, 
        " -i ", bam.file.path, " > ", 
        file.path(soup.sample.dir.paths[i], paste0(strsplit(bam.file, ".bam")[[1]], ".readistr.txt\n")), sep="")

    ### print rRNA % script
    cmd.out <- paste(cmd.out, samtools.path, " view -@ ", num.cores, " -hb -L ", hgrna.bedfile.path, " ", bam.file.path,
        "| ", samtools.path, " view -@ ", num.cores, " -hb -c > ", 
        file.path(soup.sample.dir.paths[i], paste0(strsplit(bam.file, ".bam")[[1]], ".rrnacont.txt\n")), sep="")

    cat(cmd.out,file=script.filepath,append=F)
    cat("sbatch ",script.filepath,"\n",sep="",file=jobsub.filepath, append=TRUE)
}
system(paste("chmod 700 ",jobsub.filepath,sep=""))