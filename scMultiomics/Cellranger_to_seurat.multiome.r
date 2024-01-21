###############################################
#######################################
########################## Post-alignment processing/filtering of snMultiome data
###################### Script takes cellranger output from "cellranger-arc count" function
###################### and summarizes RNA/ATAC into single seurat object
###################### Input: argument 1: path to generic folder containing cellranger output for each sample
###################### Input: argument 2: path to output folder of the results - creates analysis subfolder in output folder
###################### Output: list with filtered seurat objects
srun --mem=32GB --time=24:00:00 --pty --cpus-per-task=4 --x11 -p pi_kaminski,interactive bash
module load miniconda/4.12.0
# source activate /gpfs/gibbs/pi/kaminski/ab3478/conda_envs/agshin_R4env
source activate agshin_R4env

R

library(data.table)
library(utils)
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)
library(scales)
# library(future)
# plan("multicore", workers=4)
# options(future.globals.maxSize= 8000 * 1024^2)

args <- commandArgs(trailingOnly=TRUE)
args[1] <- "/vast/palmer/scratch/kaminski/ab3478/snMultiome_analysis/Cellranger/sample_out"
args[2] <- "/vast/palmer/scratch/kaminski/ab3478/snMultiome_analysis/Cellranger/"

### set argument 1 as top level directory
top.level.dir <- args[1]
if (dir.exists(top.level.dir) == FALSE){
    stop("No cellranger output was found\n", call = FALSE)
}


### list directories with cellranger output for each sample
sample.dirs <- list.dirs(top.level.dir, recursive=FALSE)
### list IDs of each sample
sample.names <- basename(sample.dirs)

### get gene annotations for hg38 to add gene info for the peaks
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
annotation <- renameSeqlevels(annotation, mapSeqlevels(seqlevels(annotation), "UCSC"))
genome(annotation) <- "hg38"


### Generate seurat objects for each sample
out.seurat.multiome <- list()
for (i in 1:length(sample.names)){
    ### Record start time to create combined Multiome Seurat
    start.time <- Sys.time()

    ### Set parent directory for RNA and ATAC
    parent.directory <- sample.dirs[i]
    setwd(parent.directory)

    ### Extract all unfiltered data from multiome cellranger output
    counts <- Read10X_h5(file.path(parent.directory, "outs", "raw_feature_bc_matrix.h5"))

    ### Add fragment file path
    fragpath <- file.path(parent.directory, "outs", "atac_fragments.tsv.gz")

    ### Generate Seurat object
    seurat.multiome <- CreateSeuratObject(
        project = sample.names[i],
        counts = counts$`Gene Expression`,
        assay = "RNA")
    seurat.multiome

    ### Add mitochondrial %
    seurat.multiome[["percent.mito"]] <- PercentageFeatureSet(seurat.multiome, pattern = "^MT-")

    ### Replace NaN to 0
    seurat.multiome[["percent.mito"]][is.na(seurat.multiome[["percent.mito"]])] <- 0

    ### Convert genome coordinates to matching format as in annotation
    ranges <- StringToGRanges(regions = rownames(x = counts$Peaks), sep = c(":", "-"))

    ### Create chromatin assay using cellranger peak calling
    chrom.assay <- CreateChromatinAssay(counts = counts$Peaks, 
        sep = c(":", "-"),
        fragments = fragpath, 
        ranges = ranges,
        min.cells = 0,
        min.features = 0,
        max.cells = NULL,
        motifs = NULL,
        annotation = annotation,
        bias = NULL,
        positionEnrichment = NULL)


    ### Subset barcodes present in fragment file from gene expression data
    seurat.multiome <- subset(seurat.multiome, cells = colnames(chrom.assay))

    ### Check if all barcodes match
    if (all(colnames(seurat.multiome) == colnames(chrom.assay))){
        message("All barcodes match\n")
    } else {
        message("Mismatched barcodes detected\n")
    }

    ### Add atac information to seurat object
    seurat.multiome[["ATAC"]] <- chrom.assay

    ### Change default to ATAC assay
    DefaultAssay(seurat.multiome) <- "ATAC"

    ### Create chromatin assay using MACS2 peak calling
    peaks <- CallPeaks(seurat.multiome, macs2.path = "/gpfs/gibbs/project/kaminski/ab3478/conda_envs/agshin_R4env/bin/macs2")

    ### Remove peaks on nonstandard chromosomes and in genomic blacklist regions
    peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
    n.chrom.peaks <- length(peaks)
    message("Number of peaks on standard chromosomes = ", n.chrom.peaks)

    ### Keep peaks in genomic blacklist regions to quantify number of fragments per cell
    blacklisted.peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = FALSE)
    n.blacklisted.peaks <- length(blacklisted.peaks)
    message("Number of peaks in blacklisted regions = ", n.blacklisted.peaks)

    ### Remove peaks in genomic blacklist regions
    peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

    ### Quantify counts in each peak located in the blacklisted region
    macs2_blacklisted_counts <- FeatureMatrix(fragments = Fragments(seurat.multiome), features = blacklisted.peaks, cells = colnames(seurat.multiome))

    ### Quantify counts in each peak
    macs2_counts <- FeatureMatrix(fragments = Fragments(seurat.multiome), 
        features = peaks, 
        cells = colnames(seurat.multiome))

    ### Create chromatin assay from MACS2 to the assay object
    peaks.assay <- CreateChromatinAssay(counts = macs2_counts, 
        fragments = fragpath,
        min.cells = 0,
        min.features = 0,
        max.cells = NULL,
        annotation = annotation)

    ### Subset cells only present in MACS2 counts assay (cells with 0 peaks are removed)
    seurat.multiome <- subset(seurat.multiome, cells = colnames(peaks.assay))

    ### Add peaks info from MACS2 to the assay object
    seurat.multiome[["peaks"]] <- peaks.assay

    ### Add number of blacklisted fragments and ratio in meta data
    macs2_blacklisted_counts <- macs2_blacklisted_counts[, colnames(peaks.assay)]

    seurat.multiome@meta.data$nCount_blacklist <- colSums(macs2_blacklisted_counts)

    seurat.multiome@meta.data$blacklist_ratio <- seurat.multiome@meta.data$nCount_blacklist/seurat.multiome@meta.data$nCount_peaks
    
    ### Rename cells: remove -1 and add sample ID
    seurat.multiome <- RenameCells(
        object = seurat.multiome,
        new.names = paste(sample.names[i], Cells(x = seurat.multiome), sep="__"))

    seurat.multiome <- RenameCells(
        object = seurat.multiome,
        new.names = gsub("-1", "", Cells(x = seurat.multiome)))

    ### identify number of peaks in each sample
    print(paste0("Number of peaks in sample ", sample.names[i], "=", dim(seurat.multiome@assays$peaks@counts)[1]))

    ### identify number of fragments overlapping peaks in each sample
    print(paste0("Number of fragments overlapping peaks in sample ", sample.names[i], "=", sum(seurat.multiome@assays$peaks@counts)))

    ### identify number of fragments overlaapings peaks in blacklisted regions
    print(paste0("Number of fragments overlapping peaks in blacklisted regions ", sample.names[i], "=", sum(seurat.multiome@meta.data$nCount_blacklist)))
    
    ### number of cells with > 500 fragments
    temp <- subset(seurat.multiome, nCount_peaks > 500)
    print(paste0("Number of cells with > 500 fragments in sample ", sample.names[i], "=", length(colnames(temp))))

    ### add to a list
    out.seurat.multiome[[i]] <- seurat.multiome

    ### record end time to create combined Multiome Seurat for single sample
    end.time <- Sys.time()

    ### time taken
    time.taken <- round(end.time - start.time,3)
    message(paste0("Time elapsed to process sample ", sample.names[i], " = ", time.taken, " minutes\n\n\n\n"))
}

# Number of peaks on standard chromosomes = 24574
# Number of peaks in blacklisted regions = 153
# Extracting reads overlapping genomic regions
#   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=05s  
# Extracting reads overlapping genomic regions
#   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=02m 37s
# Computing hash
# Checking for 177579 cell barcodes
# [1] "Number of peaks in sample 204_26=24421"
# [1] "Number of fragments overlapping peaks in sample 204_26=9952710"
# [1] "Number of fragments overlapping peaks in blacklisted regions 204_26=849030"
# [1] "Number of cells with > 500 fragments in sample 204_26=2996"
# Time elapsed to process sample 204_26 = 34.029 minutes

# Number of peaks on standard chromosomes = 69379
# Number of peaks in blacklisted regions = 164
# Extracting reads overlapping genomic regions
#   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=04s  
# Extracting reads overlapping genomic regions
#   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=06m 22s
# Computing hash
# Checking for 284455 cell barcodes
# [1] "Number of peaks in sample 222_68=69215"
# [1] "Number of fragments overlapping peaks in sample 222_68=29940775"
# [1] "Number of fragments overlapping peaks in blacklisted regions 222_68=560912"
# [1] "Number of cells with > 500 fragments in sample 222_68=5501"
# Time elapsed to process sample 222_68 = 33.057 minutes

# Number of peaks on standard chromosomes = 32870
# Number of peaks in blacklisted regions = 180
# Extracting reads overlapping genomic regions
#   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=04s  
# Extracting reads overlapping genomic regions
#   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=02m 45s
# Computing hash
# Checking for 203543 cell barcodes
# [1] "Number of peaks in sample 319_6=32690"
# [1] "Number of fragments overlapping peaks in sample 319_6=11416459"
# [1] "Number of fragments overlapping peaks in blacklisted regions 319_6=506017"
# [1] "Number of cells with > 500 fragments in sample 319_6=3804"
# Time elapsed to process sample 319_6 = 28.737 minutes

# Number of peaks on standard chromosomes = 69126
# Number of peaks in blacklisted regions = 189
# Extracting reads overlapping genomic regions
#   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=05s  
# Extracting reads overlapping genomic regions
#   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=07m 48s
# Computing hash
# Checking for 362245 cell barcodes
# [1] "Number of peaks in sample 325_70=68937"
# [1] "Number of fragments overlapping peaks in sample 325_70=42861357"
# [1] "Number of fragments overlapping peaks in blacklisted regions 325_70=553191"
# [1] "Number of cells with > 500 fragments in sample 325_70=10140"
# Time elapsed to process sample 325_70 = 39.093 minutes

### save the list with seurat multiome objects in Rds format
saveRDS(out.seurat.multiome, file = file.path(top.level.dir, "Seurat.combined.snATAC.Rds"))

### Create directory to store diagnostic plots
analysis.dir <- file.path(args[2], "analysis")
if (dir.exists(analysis.dir) == FALSE){
    dir.create(analysis.dir)
}

setwd(analysis.dir)

############################################
##########################
################ Quality Control (QC)
############ snATAC-seq
### Look at nCount_RNA, nCount_peaks, nucleosome_signal, TSS.enrichment
for (i in 1:length(out.seurat.multiome)){

    DefaultAssay(out.seurat.multiome[[i]]) <- "peaks"
    ### Estimate nucleosome signal
    out.seurat.multiome[[i]] <- NucleosomeSignal(object = out.seurat.multiome[[i]])
    ### Estimate tss enrichment
    out.seurat.multiome[[i]] <- TSSEnrichment(object = out.seurat.multiome[[i]], fast=FALSE)

    sample.name <- levels(out.seurat.multiome[[i]]@meta.data$orig.ident)
    ### Estimate high and low tss enrichment (threshold = 2)
    out.seurat.multiome[[i]]$high.tss <- ifelse(out.seurat.multiome[[i]]$TSS.enrichment > 2, 'High', 'Low')
    ### Estimate high and low nucleosome signal (threshold = 4)
    out.seurat.multiome[[i]]$nucleosome_group <- ifelse(out.seurat.multiome[[i]]$nucleosome_signal > 4, 'NS > 4', 'NS < 4')

    ### Make TSSplot
    png(paste0("TSSPlot.snmultiome.", sample.name, ".01302023.png"), res=200, unit="in", height=8, width=11)
    print(TSSPlot(out.seurat.multiome[[i]], group.by = 'high.tss') + NoLegend())
    dev.off()

    ### Make fragment histogram
    png(paste0("FragmentHistogram.snmultiome.", sample.name, ".01302023.png"), res=200, unit="in", height=8, width=11)
    print(FragmentHistogram(object = out.seurat.multiome[[i]], group.by = 'nucleosome_group'))
    dev.off()

}

### Create temporary seurat object to make violin plot
out.counts <- list()
out.metadata <- list()
for (i in 1:length(out.seurat.multiome)){
    out.counts[[i]] <- out.seurat.multiome[[i]]@assays$RNA@counts
    out.metadata[[i]] <- out.seurat.multiome[[i]]@meta.data
}

counts <- do.call(cbind, out.counts)
metadata <- do.call(rbind, out.metadata)

temp.seurat <- CreateSeuratObject(counts, meta.data = metadata)
temp.seurat
# An object of class Seurat 
# 36601 features across 1027822 samples within 1 assay 
# Active assay: RNA (36601 features, 0 variable features)

### Make violin plot
png(paste0("VlnPlot.snmultiome.snatac.01302023.png"), res=200, unit="in", height=8, width=11)
VlnPlot(object = temp.seurat, 
    features = c("nCount_RNA", "nCount_peaks", "TSS.enrichment", "nucleosome_signal", "blacklist_ratio"),
    pt.size = 0) + scale_y_continuous(limits = c(0, 2000))
dev.off()

VlnPlot(object = temp.seurat, features = c("nCount_peaks"), pt.size = 0) + scale_y_continuous(limits = c(0, 1000))

hist(out.seurat.multiome[[1]]@meta.data$nCount_peaks)

########################## knee plot function
betterBarcode_rank_plot <- function(sce, nBarcodes, sample){
    counts <- sce@meta.data$nCount_RNA[order(sce@meta.data$nCount_RNA, decreasing=TRUE)]
    counts <- sort(counts, decreasing = TRUE)
    ranks <- seq(length(counts))
    datf <- data.frame(Rank = ranks, Count = counts)
    datf <- datf[!is.na(datf[, 2]), , drop = FALSE]

    counts[duplicated(counts)] <- NA
    ranks[duplicated(counts)] <- NA
    datf.plot <- data.frame(Rank = ranks, Count = counts)
    datf.plot <- datf.plot[!is.na(datf.plot[, 2]), , drop = FALSE]
    png(paste("Elbow.nUMI.total.", sample, ".png", sep=""), res=200, unit="in", height=8, width=11)
    print({
        p <- ggplot(datf.plot, aes_string(x = "Rank", y = "Count")) +
        geom_point() +
        scale_x_continuous(trans = "log10", labels = comma,
            breaks=c(10, 100, 500, 1000, 5000, 10000, 30000, 50000, 100000, 500000),
            limits=c(1, 500000)) +
        scale_y_continuous(name = "Droplet size", trans = "log10",labels = comma,
            breaks=c(1, 10, 100, 150, 200, 500, 1000, 5000, 10000, 25000, 50000),
            limits=c(1, 50000)) +
            geom_vline(xintercept=nBarcodes, linetype="dashed", color="red")  +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle(sample)
    })
    dev.off()
    return(datf)
}

### Create directory to store output processed seurat objects
processed_sample.dir <- file.path(args[2], "processed_data")
if (dir.exists(processed_sample.dir) == FALSE){
    dir.create(processed_sample.dir)
}

setwd(processed_sample.dir)

#### use knee plot function to extract 30000 barcodes with highest nCount_RNA
nBarcodes <- 30000
mega.seurat.list <- list()
for (i in 1:length(out.seurat.multiome)){
    temp.seurat <- out.seurat.multiome[[i]]
    rank.matrix <- betterBarcode_rank_plot(temp.seurat, nBarcodes, unique(temp.seurat@meta.data$orig.ident))
    umi.threshold <- rank.matrix[rank.matrix$Rank == nBarcodes,]$Count
    temp.seurat <- subset(temp.seurat, nCount_RNA >= umi.threshold)
    mega.seurat.list[[i]] <- temp.seurat
    names(mega.seurat.list)[i] <- levels(mega.seurat.list[[i]])
    print(mega.seurat.list[[i]])

    ### number of cells with > 500 fragments
    temp <- subset(temp.seurat, nCount_peaks > 500)
    print(paste0("Number of cells with > 500 fragments in sample ", sample.names[i], "=", length(colnames(temp))))

}

# An object of class Seurat 
# 95422 features across 30181 samples within 3 assays 
# Active assay: peaks (24421 features, 0 variable features)
#  2 other assays present: RNA, ATAC
# [1] "Number of cells with > 500 fragments in sample 204_26=2993"

# An object of class Seurat 
# 165150 features across 31046 samples within 3 assays 
# Active assay: peaks (69215 features, 0 variable features)
#  2 other assays present: RNA, ATAC
# [1] "Number of cells with > 500 fragments in sample 222_68=5334"

# An object of class Seurat 
# 103077 features across 30006 samples within 3 assays 
# Active assay: peaks (32690 features, 0 variable features)
#  2 other assays present: RNA, ATAC
# [1] "Number of cells with > 500 fragments in sample 319_6=3774"

# An object of class Seurat 
# 171391 features across 30248 samples within 3 assays 
# Active assay: peaks (68937 features, 0 variable features)
#  2 other assays present: RNA, ATAC
# [1] "Number of cells with > 500 fragments in sample 325_70=9968"

save(mega.seurat.list, file = "Seurat.multiome.013123.Robj")

















