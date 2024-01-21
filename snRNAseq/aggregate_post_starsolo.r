#aggregating data after STARSolo


#  saved as aggregate_STARSolo.R

### request small amount of resources
srun --mem=32GB --time=6:00:00 --pty --cpus-per-task=4 --x11 -p bigmem,interactive bash

## load conda to enter R environment
library(Matrix)
library(tidyverse)
library(ggplot2)
library(scales)
library(cowplot)
library(Seurat)
library(RColorBrewer)
library(viridis)
library(viridisLite)
library(circlize)
options(stringsAsFactor=FALSE)


########################## knee plot function
betterBarcode_rank_plot <- function(sce, nBarcodes, sample){
    counts <- sce$metadata$nUMI[order(sce$metadata$nUMI, decreasing=TRUE)]
    counts <- sort(counts, decreasing = TRUE)
    ranks <- seq(length(counts))
    datf <- data.frame(Rank = ranks, Count = counts)
    datf <- datf[!is.na(datf[, 2]), , drop = FALSE]

    counts[duplicated(counts)] <- NA
    ranks[duplicated(counts)] <- NA
    datf.plot <- data.frame(Rank = ranks, Count = counts)
    datf.plot <- datf.plot[!is.na(datf.plot[, 2]), , drop = FALSE]
    png(paste("Elbow.nUMI.total.", sce$metadata$library.ident, ".png", sep=""), res=200, unit="in", height=8, width=11)
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

sce <- output
sample <- output$metadata$library.ident
nBarcodes <- 30000
## this the parent directory for the sample outputs
topLevel.dir <- "/vast/palmer/home.grace/aj597/scratch60/Ed_cutadap/sample_out"
### we will work in analysis folder
setwd("/vast/palmer/home.grace/aj597/scratch60")

sample.paths <- list.files(topLevel.dir, full.names=TRUE)
sample.names <- basename(sample.paths)

##### loop over samples, pull the filtered matrix, append sample identifier, collapse the matrices
### make an empty list



list.counts <- list()
list.metadata <- list()

collectSTARsolo <- function(sample.path, nBarcodes){
    sample.name <- basename(sample.path)
    cat("collecting \'GeneFull\' count information\n")
#### get the gene-full information first
    geneFull.path <- file.path(sample.path, "Solo.out", "GeneFull","raw")
    gf.temp.barcodes <- read.table(file.path(geneFull.path, "barcodes.tsv"), sep="\t", header=FALSE)[,1]
    gf.temp.genes <- as.vector(read.table(file.path(geneFull.path, "features.tsv"), sep="\t", header=FALSE)$V2)
    gf.counts <-  Matrix::readMM(file.path(geneFull.path, "matrix.mtx"))
### put them together
    rownames(gf.counts) <- gf.temp.genes
    colnames(gf.counts) <- gf.temp.barcodes
############################################
########### gathering the rest of the stuff
    cat("collecting \'Gene\' and \'Velocyto\' splice information\n")
    gene.path <- file.path(sample.path, "Solo.out", "Gene","raw")
    g.temp.barcodes <- read.table(file.path(gene.path, "barcodes.tsv"), sep="\t", header=FALSE)[,1]
    g.temp.genes <- as.vector(read.table(file.path(gene.path, "features.tsv"), sep="\t", header=FALSE)$V2)
    g.counts <-  Matrix::readMM(file.path(gene.path, "matrix.mtx"))
### put them together
    rownames(g.counts) <- g.temp.genes
    colnames(g.counts) <- g.temp.barcodes
#### get % umi in Gene vs GeneFull
    percSpliced.g <- 100*(Matrix::colSums(g.counts)/Matrix::colSums(gf.counts))
## replace NA w/ zero
    percSpliced.g[is.nan(percSpliced.g)] <- 0
## replace Inf w/ 100
    percSpliced.g[is.infinite(percSpliced.g)] <- 100
################## get the scVelo spliced/unspliced estimate of spliced fraction
    velocyto.path <- file.path(sample.path, "Solo.out", "Velocyto","raw")
    spliced.mtx <- readr::read_delim(file.path(velocyto.path, "spliced.mtx"), col_names = FALSE, show_col_types=FALSE, delim = ' ', skip = 3)
    unspliced.mtx <- readr::read_delim(file.path(velocyto.path, "unspliced.mtx"), col_names = FALSE, show_col_types=FALSE, delim = ' ', skip = 3)
    ambiguous.mtx <- readr::read_delim(file.path(velocyto.path, "ambiguous.mtx"), col_names = FALSE, show_col_types=FALSE, delim = ' ', skip = 3)
    # velo.t <- readr::read_delim(file.path(velocyto.path, "matrix.mtx"), col_names = FALSE, delim = ' ', skip = 3)
    v.temp.barcodes <- read.table(file.path(velocyto.path, "barcodes.tsv"), sep="\t", header=FALSE)[,1]
#### split the three levels of information into separate matrices
    spliced.mtx <- Matrix::sparseMatrix(i = spliced.mtx$X1, j = spliced.mtx$X2, x = spliced.mtx$X3, dims = dim(gf.counts))
    unspliced.mtx <- Matrix::sparseMatrix(i = unspliced.mtx$X1, j = unspliced.mtx$X2, x = unspliced.mtx$X3, dims = dim(gf.counts))
    ambiguous.mtx <- Matrix::sparseMatrix(i = ambiguous.mtx$X1, j = ambiguous.mtx$X2, x = ambiguous.mtx$X3, dims = dim(gf.counts))
#### get an aggregated colsum of all barcodes
    all.colsums <- Matrix::colSums(spliced.mtx) + Matrix::colSums(unspliced.mtx) + Matrix::colSums(ambiguous.mtx)
#### assign the barcode names to the values
    percSpliced.v <- 100*(Matrix::colSums(spliced.mtx)/all.colsums)
    names(percSpliced.v) <- v.temp.barcodes
## replace NA w/ zero and Inf w/ 100, assign names
    percSpliced.v[is.nan(percSpliced.v)] <- 0
    percSpliced.v[is.infinite(percSpliced.v)] <- 100
    percUnSpliced.v <- 100*(Matrix::colSums(unspliced.mtx)/all.colsums)
    names(percUnSpliced.v) <- v.temp.barcodes
## replace NA w/ zero and Inf w/ 100, assign names
    percUnSpliced.v[is.nan(percUnSpliced.v)] <- 0
    percUnSpliced.v[is.infinite(percUnSpliced.v)] <- 100
## replace NA w/ zero and Inf w/ 100, assign names
    percAmbiguous.v <- 100*(Matrix::colSums(ambiguous.mtx)/all.colsums)
    percAmbiguous.v[is.nan(percAmbiguous.v)] <- 0
    percAmbiguous.v[is.infinite(percAmbiguous.v)] <- 100
    names(percAmbiguous.v) <- v.temp.barcodes
#### now only keep the barcodes selected from nBarcodes
#     percSpliced.v <- percSpliced.v[barcodes.keep]
#     percUnSpliced.v <- percUnSpliced.v[barcodes.keep]
#     percAmbiguous.v <- percAmbiguous.v[barcodes.keep]
#### collect the remaining metadata values of interest
    gf.nUMI <- Matrix::colSums(gf.counts)
    gf.nGenes <- Matrix::colSums(gf.counts > 0)
    mito.genes <- grep("^mt-", rownames(gf.counts), value=TRUE)
    percMito <- 100*(Matrix::colSums(gf.counts[mito.genes,])/gf.nUMI)
    percMito[is.nan(percMito)] <- 0
    percMito[is.infinite(percMito)] <- 100
    percMALAT1 <- 100*(gf.counts["Malat1",]/gf.nUMI)
    percMALAT1[is.nan(percMALAT1)] <- 0
    percMALAT1[is.infinite(percMALAT1)] <- 100
### output each iteration into spot on list
    metadata <- data.frame(barcode=colnames(gf.counts),
        nUMI=gf.nUMI[colnames(gf.counts)], nGene=gf.nGenes[colnames(gf.counts)],
        percMito=percMito[colnames(gf.counts)], percMALAT1=percMALAT1[colnames(gf.counts)],
        percSpliced.g=percSpliced.g, percSpliced.v=percSpliced.v[colnames(gf.counts)],
        percUnspliced.v=percUnSpliced.v[colnames(gf.counts)], percAmbiguous.v=percAmbiguous.v[colnames(gf.counts)],
        library.ident=rep.int(sample.name, times=length(colnames(gf.counts))))
    out.list <- list(gf.counts, metadata)
    names(out.list) <- c("counts", "metadata")
    rank.matrix <- betterBarcode_rank_plot(out.list, nBarcodes, out.list$metadata$library.ident)
    umi.threshold <- rank.matrix[rank.matrix$Rank == nBarcodes,]$Count
    cat("nUMI threshold for library ", sample.names[i] , " = ", umi.threshold, "\n")
    rownames(out.list$metadata) <- paste(out.list$metadata$library.ident, "__", rownames(out.list$metadata),sep="")
    colnames(out.list$counts) <- rownames(out.list$metadata)
    out.list$metadata=out.list$metadata[out.list$metadata$nUMI >= umi.threshold,]
    cat("Number of barcodes extracted from library ", sample.names[i] , " = ", nrow(out.list$metadata), "\n")
    out.list$counts=out.list$counts[,colnames(out.list$counts) %in% rownames(out.list$metadata)]
    return(out.list)
}


library.keep<- c("OF", "F_HOX", "M_HOX", "YF", "YM", "OM", "F_NOX", "M_NOX")
sample.names= sample.names[sample.names %in% library.keep]
sample.paths=sample.paths[basename(sample.paths) %in% library.keep]

for (i in 1:length(sample.paths)) {
  output <- collectSTARsolo(sample.paths[i], nBarcodes = 30000)
  list.counts[[i]] <- output$counts
  list.metadata[[i]] <- output$metadata
}

##Cresate big seurat Object
count.matrix=do.call(cbind, list.counts)
metadata.matrix = do.call(rbind, list.metadata)

dim(count.matrix)
# [1]  53647 244066
dim(metadata.matrix)
# [1] 244066     10
sub.ed <- CreateSeuratObject(counts=count.matrix, meta.data=metadata.matrix, names.delim="__")

####### Add batch column
batch1<-c("F_HOX", "F_NOX", "M_HOX","M_NOX")
batch2<-c("OF" ,"OM" ,"YF" ,"YM")
batch.list=c(rep(NA, dim(metadata.matrix)[1]))
for (i in c(batch1,batch2)) {
  if (i %in% batch1){
    positions <- which(sub.ed@meta.data$orig.ident == i)
    batch.list[positions] <- "batch1"
  }
  if (i %in% batch2) {
    positions <- which(sub.ed@meta.data$orig.ident == i)
    batch.list[positions] <- "batch2"
  }
}

sub.ed@meta.data$batch <- batch.list



#################
########## Make diagnostic plots
## Make log10nUMI column
sub.ed@meta.data$log10nUMI <- log10(sub.ed@meta.data$nUMI)

########## Scatter plots
###### nUMI vs percSpliced
png("ScatterPlot.nUMI.percSpliced.ED.png", res=200, unit="in", height=8, width=11)
ggplot(sub.ed@meta.data, aes(x=percSpliced.v, y=log10nUMI)) +
    geom_point(aes(color=batch)) +
    xlab("Percent Spliced") +
    ylab("nUMI") +
    scale_fill_manual(values = c("batch1"="yellow", "batch2"="red")) +
    scale_x_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100)) +
    # scale_color_viridis(discrete = TRUE, option = "C")+
    # scale_fill_viridis(discrete = TRUE) +
    theme(panel.background = element_rect(fill = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
dev.off()


###### nUMI vs percUnspliced

png("ScatterPlot.nUMI.percUnspliced.ED.png", res=200, unit="in", height=8, width=11)
ggplot(sub.ed@meta.data, aes(x=percUnspliced.v, y=log10nUMI)) +
    geom_point(aes(color=batch)) +
    xlab("Percent Intronic") +
    ylab("nUMI") +
    scale_fill_manual(values = c("batch1"="yellow", "batch2"="red")) +
    scale_x_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100)) +
    geom_vline(xintercept =7,color="white") +
    # scale_color_viridis(discrete = TRUE, option = "C")+
    # scale_fill_viridis(discrete = TRUE) +
    theme(panel.background = element_rect(fill = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
dev.off()


###### nUMI vs Mitochondrial
png("ScatterPlot.nUMI.percMito.ED.png", res=200, unit="in", height=8, width=11)
ggplot(sub.ed@meta.data, aes(x=percMito, y=log10nUMI)) +
    geom_point(aes(color=batch)) +
    xlab("Percent Mitochondrial") +
    ylab("nUMI") +
    scale_fill_manual(values = c("batch1"="yellow", "batch2"="red")) +
    # scale_color_viridis(discrete = TRUE, option = "D")+
    # scale_fill_viridis(discrete = TRUE) +
    theme(panel.background = element_rect(fill = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
dev.off()

#####nUMIvs genes
png("ScatterPlot.nUMI.gene.ED.png", res=200, unit="in", height=8, width=11)
ggplot(sub.ed@meta.data, aes(x=nGene , y=log10nUMI)) +
    geom_point(aes(color=batch)) +
    xlab("Gene") +
    ylab("nUMI") +
    scale_fill_manual(values = c("batch1"="yellow", "batch2"="red")) +
    #scale_x_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100)) +
    geom_vline(xintercept =300,color="white") +
    # scale_color_viridis(discrete = TRUE, option = "C")+
    #scale_fill_viridis(discrete = TRUE) +
    theme(panel.background = element_rect(fill = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
dev.off()


#####Caraterize data unsliced<8%
##Metadata add 1 if for barcodes with >8%
sub.ed@meta.data$statut_unsp <- c(rep(NA,dim(metadata.matrix)[1]))
sub.ed@meta.data[which( sub.ed@meta.data$percUnspliced.v <= 8),"statut_unsp"] = "0"
sub.ed@meta.data[which( sub.ed@meta.data$percUnspliced.v > 8),"statut_unsp"] = "1"
##Make the figure
png("ScatterPlot.nUMI.gene.unsp.ED.png", res=200, unit="in", height=8, width=11)
ggplot(sub.ed@meta.data, aes(x=nGene , y=log10nUMI)) +
    geom_point(aes(color=statut_unsp)) +
    xlab("Gene") +
    ylab("nUMI") +
    scale_fill_manual(values = c("0"="yellow", "1"="red")) +
    #scale_x_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100)) +
    geom_vline(xintercept =300,color="white") +
    # scale_color_viridis(discrete = TRUE, option = "C")+
    #scale_fill_viridis(discrete = TRUE) +
    theme(panel.background = element_rect(fill = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
dev.off()

png("ScatterPlot.nUMI.percMito.unsp.ED.png", res=200, unit="in", height=8, width=11)
ggplot(sub.ed@meta.data, aes(x=percMito, y=log10nUMI)) +
    geom_point(aes(color=statut_unsp)) +
    xlab("Percent Mitochondrial") +
    ylab("nUMI") +
    scale_fill_manual(values = c("0"="yellow", "1"="red")) +
    # scale_color_viridis(discrete = TRUE, option = "D")+
    # scale_fill_viridis(discrete = TRUE) +
    theme(panel.background = element_rect(fill = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
dev.off()
