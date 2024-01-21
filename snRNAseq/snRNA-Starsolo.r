## argument 1: path to parent directory for the project
## eg: Rscript 10Xv3.STARsolo.snLeuven.R /home/fas/kaminski/tsa27/scratch60/snRNA_leuven
#####################################################################################################################################################
########## This Rscript is designed to use as input: the same directory labeled as "output.parent.dir" in the cutadapt Rscript
########## which contains a directory titles "sample_out" and "cutadapt_scripts"
########## The cutadapt job must be run beforehand to generate "trimmed.R1.fastq" and "trimmed.R2.fastq" in the
########## "trimmed_merged_fastq" subdirectory in each sample's output directory
##########
########## The purpose of the script is to generate a sample-specific .sh rocessing each sample through STARsolo
##########
###### Summary:
###### 1 Identify all sample folders
#loop{ 2 generate script for STAR solo}
#####################################################################################################################################################
# 8 threads
# 60GB memory
# for V2 and V3 you need to change the length of the UMI and the whitelist

args <- commandArgs(trailingOnly=TRUE)
args[1] <- "/vast/palmer/scratch/kaminski/ab3478/snRNA_Multiome"
# args[1] is input parent directory

sample.parent.dir <- file.path(args[1], "sample_out")
scripts.parent.dir <- file.path(args[1], "scripts")
STARsolo.scripts.dir <- file.path(scripts.parent.dir, "STARsolo")
if(file.exists(STARsolo.scripts.dir)==F){
            dir.create(STARsolo.scripts.dir)
}

soup.sample.names <- list.files(file.path(sample.parent.dir))
soup.sample.dir.paths <- file.path(sample.parent.dir, soup.sample.names)
jobsub.filepath <- file.path(STARsolo.scripts.dir,"STARsolo.jobsub.bat")

STAR.path <- "/gpfs/gibbs/pi/kaminski/public/softwares/STAR-2.7.6a/bin/Linux_x86_64_static/STAR"
genome.dir <- "/gpfs/gibbs/pi/kaminski/public/Backup/Jonas/genome_index_human/GENCODE_release37_GRCh38.p13/STARindex_noGTF/"
GTF.path <-  "/gpfs/gibbs/pi/kaminski/public/Backup/Jonas/genome_index_human/GENCODE_release37_GRCh38.p13/GTF/gencode.v37.annotation.gtf"

## for Single Cell Multiome (ATAC+GEX) v1
CBwhitelist.path <- "/vast/palmer/scratch/kaminski/ab3478/10x_White_Lists/737K-jag-v1.txt"


STAR.options <- paste("--soloType CB_UMI_Simple --readFilesCommand zcat --soloFeatures Gene GeneFull SJ Velocyto --runThreadN 8 ",
    "--soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 --twopassMode Basic --soloStrand Forward", sep="") # --soloUMIlen 10 for V2, --soloUMIlen 12 for V3

## Wipe the old jobsub.bat file
cat("",sep="",file=jobsub.filepath, append=FALSE)
## Loop over samples to extract trimmed R1 and R2 1:length(soup.sample.dir.paths)
for(i in 1:length(soup.sample.names)){
      trimmed.folder.paths <- file.path(soup.sample.dir.paths[i], "trimmed_merged_fastq")
#### Distinguish Read1 and Read2
    read1.fastq.path <- file.path(trimmed.folder.paths , "trimmed.R1.fastq.gz")
    read2.fastq.path <- file.path(trimmed.folder.paths , "trimmed.R2.fastq.gz")
    script.filepath <- file.path(STARsolo.scripts.dir, paste(soup.sample.names[i], "_STARsolo.sh", sep = ""))
##### print the script ##  -N ",num.reads.per.cell,";
    cmd.out <- NULL
    cmd.out <- paste("#!/bin/bash\n")
    cmd.out <- paste(cmd.out,"#SBATCH --time=6:00:00 -p day,pi_kaminski --ntasks=1 --cpus-per-task=8",
        " --mem=60000M --job-name=", paste0(soup.sample.names[i], "_STARsolo"),
        " -o ",script.filepath,".o%J -e ",script.filepath,".e%J\n",sep="")
    cmd.out <- paste(cmd.out, STAR.path, " ", STAR.options, " --soloCBwhitelist ", CBwhitelist.path,
        " --genomeDir ", genome.dir, " --sjdbGTFfile ", GTF.path,
        " --readFilesIn ", read2.fastq.path, " ", read1.fastq.path,
        " --outFileNamePrefix ", soup.sample.dir.paths[i], "/", sep="")
    cat(cmd.out,file=script.filepath,append=F)
    cat("sbatch ",script.filepath,"\n",sep="",file=jobsub.filepath, append=TRUE)
}
system(paste("chmod 700 ",jobsub.filepath,sep=""))