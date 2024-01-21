#!/bin/bash

############## Generating RNA-Seq count matrix
############## This file inputs bulk RNA-seq BAM files and generates counts matrix output using featureCounts software
############## Input: RNA-Seq BAM file
############## Output: Analysis-ready RNA COUNT file and the SUMMARY


### Number of available cores
num.cores <- 8

## Set input directory with BAM files
input.parent.dir <- "/home/abalay/scratch/PD1I_datasets/processed/sample_out"

## Set output directory
output.parent.dir <- "/home/abalay/scratch/PD1I_datasets/Transcriptomics"
if (dir.exists(output.parent.dir) == FALSE){
    message("Creating output directory ...", "\n")
    dir.create(output.parent.dir)
}

### Create directory to store counts matrices
counts.dir <- file.path(output.parent.dir, "Counts")
if (dir.exists(counts.dir) == FALSE){
    message("Creating counts directory ...", "\n")
    dir.create(counts.dir)
}

## Load genome annotation file (in gtf format)
gtf.path <- "/home/abalay/scratch/Pipelines_and_Files/Files_for_scripts/gencode.v37.annotation.gtf"

## Load sample metadata 
metadata.path <- file.path("/home/abalay/scratch/PD1I_datasets", "PD1I_metadata.txt")
pd1l.metadata <- read.table(metadata.path, sep = "\t", header=TRUE)

## Create count matrices for different library designs (single/paired and stranded/unstranded)

### label samples with high rRNA read % (i.e. 20%)
excluded.samples <- c("SRR21531011", "SRR3184285", "SRR3184288", "SRR3184298", "SRR5088856", "SRR5088872")

soup.sample.names <- list.files(input.parent.dir)
soup.sample.dir.paths <- file.path(input.parent.dir, soup.sample.names)

### append sample names to vectors based on library design
library.list <- list()

library.list[[1]] = library.list[[2]] = library.list[[3]] = library.list[[4]] = c("")

names(library.list) <- c("ss", "su", "ps", "pu")

for (i in 1:length(soup.sample.names)){
### set strandedness and library layout of sample
    strand.info <- pd1l.metadata$stranded.unstranded_rnaseq[pd1l.metadata$Run == soup.sample.names[i]]
    lib.layout.info <- pd1l.metadata$LibraryLayout[pd1l.metadata$Run == soup.sample.names[i]]

### set path to bam file
    bam.path <- list.files(soup.sample.dir.paths[i], pattern = ".bam", full.names=TRUE)[1]
    if (lib.layout.info == "SINGLE" & strand.info == "stranded"){
        library.list$ss <- append(library.list$ss, bam.path)
    }
    if (lib.layout.info == "SINGLE" & strand.info == "unstranded"){
        library.list$su <- append(library.list$su, bam.path)
    }
    if (lib.layout.info == "PAIRED" & strand.info == "stranded"){
        library.list$ps <- append(library.list$ps, bam.path)
    }
    if (lib.layout.info == "PAIRED" & strand.info == "unstranded"){
        library.list$pu <- append(library.list$pu, bam.path)
    }
}

### generate count matrices
for (j in 1:length(library.list)){
    if (length(library.list[[j]]) > 1){
        if (names(library.list)[[j]] == "ss"){
            ss.samples <- paste(library.list[[j]], collapse = " ")
            system(paste("featureCounts", "-T", num.cores, "-s 2", "-t exon", "-g gene_name", "-a", gtf.path, "-o",
                file.path(counts.dir, "SS.counts.txt"),
                ss.samples, sep=" "))
        }
        if (names(library.list)[[j]] == "su"){
            su.samples <- paste(library.list[[j]], collapse = " ")
            system(paste("featureCounts", "-T", num.cores, "-s 0", "-t exon", "-g gene_name", "-a", gtf.path, "-o", 
                file.path(counts.dir, "SU.counts.txt"),
                su.samples, sep=" "))
        }
        if (names(library.list)[[j]] == "ps"){
            ps.samples <- paste(library.list[[j]], collapse = " ")
            system(paste("featureCounts", "-p --countReadPairs", "-T", num.cores, "-s 2", "-t exon", "-g gene_name", "-a", gtf.path, 
                "-o", file.path(counts.dir, "PS.counts.txt"), ps.samples, sep=" "))
        }
        if (names(library.list)[[j]] == "pu"){
            pu.samples <- paste(library.list[[j]], collapse = " ")
            system(paste("featureCounts", "-p --countReadPairs", "-T", num.cores, "-s 0", "-t exon", "-g gene_name", "-a", gtf.path, 
                "-o", file.path(counts.dir, "PU.counts.txt"), pu.samples, sep=" "))
        }
    }
}

### list all counts matrices
counts.file.paths <- list.files(counts.dir, pattern = "\\.counts.txt$", full.names=TRUE)

list.counts <- list()
### merge all counts matrices into single count matrix
for (i in 1:length(counts.file.paths)){
    counts <- read.table(counts.file.paths[i], sep = "\t", header = TRUE, row.names = "Geneid")

    columns.delete <- c("Chr", "Start", "End", "Strand", "Length")
    counts <- counts[, !colnames(counts) %in% columns.delete]

    list.counts[[i]] <- counts
}

big.count.matrix <- do.call(cbind, list.counts)
dim(big.count.matrix)
# [1] 59409   175

### save merged count matrix
write.csv(big.count.matrix, file = file.path(counts.dir, "big.featureCounts.matrix.csv"), row.names = TRUE)












