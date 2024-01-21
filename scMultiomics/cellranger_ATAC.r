## argument 1: directory with raw single-nuclei RNA-SEQ fastq files
## argument 2: directory with raw single-nuclei ATAC-SEQ fastq files
## argument 3: output directory
## eg: Rscript cellranger_ATAC.r -p /gpfs/gibbs/pi/kaminski/public/Backup/Agshin/snMultiome -o /gpfs/gibbs/pi/kaminski/public/Backup/Agshin/snMultiome_analysis
#####################################################################################################################################################
########## The purpose of the script is to generate a sample-specific shell file processing each sample through Cellranger
##########
###### Summary: 
###### 1 Identify all sample folders
#loop{ 2 generate script for Cellranger}
#####################################################################################################################################################
# 8 threads
# 100GB memory

## Number of cores to multithread
num.cores <- 8

### install required libraries
require("optparse")

### set the arguments
option_list = list(
    make_option(c("-r", "--rnadir"), type="character", default=NULL, 
              help="directory with raw snRNA-seq fastq files", metavar="character"),
    make_option(c("-a", "--atacdir"), type="character", default=NULL, 
              help="directory with raw snATAC-seq fastq files", metavar="character"),
    make_option(c("-o", "--output-dir"), type="character", default=NULL, 
              help="output directory for Cellranger output", metavar="character")); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

### Return error and cancel immediately if input directory doesn't exist
input.parent.dir <- opt$p
if(file.exists(input.parent.dir)==F){
    cat("Error: Input parent directory does not exist... Try again.\n")
    exit()
}
# specify path to RNA and ATAC batches
raw.RNA.file.dir <- opt$'rnadir'
raw.ATAC.file.dir <- opt$'atacdir'


## Generate the output directory structure
output.parent.dir <- opt$o
if(file.exists(output.parent.dir)==F){
    cat("Generating top level directory for output...\n")
    dir.create(output.parent.dir)
}

# create upper level directory for output
top.level.dir <- file.path(output.parent.dir, "Cellranger")
if(file.exists(top.level.dir)==F){
    dir.create(top.level.dir)
}

# create subfolder for sample output
sample.parent.dir <- file.path(top.level.dir, "sample_out")
if(file.exists(sample.parent.dir)==F){
    dir.create(sample.parent.dir)
}

# create subfolder for script output
scripts.parent.dir <- file.path(top.level.dir, "scripts")
if(file.exists(scripts.parent.dir)==F){
    dir.create(scripts.parent.dir)
}

# create directory to output cellranger scripts
cellranger.scripts.dir <- file.path(scripts.parent.dir, "cellranger-arc-2.0.2")
if(file.exists(cellranger.scripts.dir)==F){
    dir.create(cellranger.scripts.dir)
}

# batch level information 
soup.ATAC.batch.names <- list.files(raw.ATAC.file.dir)
soup.ATAC.batch.dir.paths <- file.path(raw.ATAC.file.dir, soup.ATAC.batch.names)

# create file with whole cellranger script
jobsub.filepath <- file.path(cellranger.scripts.dir,"cellranger-arc.jobsub.bat")

# add cellranger path
cellranger.path <- "/vast/palmer/scratch/kaminski/ab3478/softwares/cellranger-arc-2.0.2/cellranger-arc"
# add genome path
genome.dir <- "/vast/palmer/scratch/kaminski/ab3478/Refgenome/refdata-cellranger-arc-GRCh38-2020-A-2.0.0"
# add 10x whitelist path
CBwhitelist.path <- "/vast/palmer/scratch/kaminski/ab3478/softwares/cellranger-arc-2.0.2/lib/python/atac/barcodes/737K-arc-v1.txt"

## Wipe the old jobsub.bat file
cat("",sep="",file=jobsub.filepath, append=FALSE)
## Loop over samples to run cellranger in multiome mode
for(i in 1:length(soup.ATAC.batch.names)){
    soup.ATAC.sample.names <- list.files(soup.ATAC.batch.dir.paths[i])
    soup.ATAC.sample.dir.paths <- file.path(soup.ATAC.batch.dir.paths[i],soup.ATAC.sample.names)
    for (j in 1:length(soup.ATAC.sample.names)){
        sample.output.dir <- file.path(sample.parent.dir, soup.ATAC.sample.names[j])
        if(file.exists(sample.output.dir)==F){
            dir.create(sample.output.dir)
        }
        script.filepath <- file.path(cellranger.scripts.dir, paste0(soup.ATAC.sample.names[j], "_cellranger.sh"))
        library.file <- file.path(input.parent.dir, paste0("library_", soup.ATAC.sample.names[j], ".csv"))
        ##### print the script ##  -N ",num.reads.per.cell,";
        cmd.out <- NULL
        cmd.out <- paste("#!/bin/bash\n")
        cmd.out <- paste(cmd.out,"#SBATCH --time=48:00:00 -p week,pi_kaminski --ntasks=1 --cpus-per-task=8",
            " --mem=100G --job-name=", paste0(soup.ATAC.sample.names[j], "_cellranger"),
            " -o ",script.filepath,".o%J -e ",script.filepath,".e%J\n",sep="")
        cmd.out <- paste(cmd.out, "cd ", sample.output.dir, "\n", sep="")
        cmd.out <- paste(cmd.out, cellranger.path, " count ",
            "--id=",soup.ATAC.sample.names[j],
            " --reference=", genome.dir,
            " --libraries=", library.file,
            " --localcores=8 ",
            "--localmem=60", "\n\n", sep="")  
        cat(cmd.out,file=script.filepath,append=F)
        cat("sbatch ",script.filepath,"\n",sep="",file=jobsub.filepath, append=TRUE)
    }
}

system(paste("chmod 700 ",jobsub.filepath,sep=""))