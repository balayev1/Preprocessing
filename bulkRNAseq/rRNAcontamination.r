#############  Mapping Quality Statistics Analysis
##### rRNA contamination level RNA-Seq

library(ggplot2)
#### Upload stats and number of rRNA-aligned reads text files

### Working directory
wd <- "/home/abalay/scratch/PD1I_datasets/processed"

### List all folders
dirs <- list.dirs(file.path(wd, "sample_out"), recursive = FALSE)

### Select only sample folders
# removed.folders <- c("/home/cluster/abalay/scratch/Cell_Lines/Scripts", "/home/cluster/abalay/scratch/Cell_Lines/multiqc_reports")
# dirs <- dirs[!dirs %in% removed.folders]


### Add sample name, rRNA read number and total mapped read number
sample_name = c(); rrna = c(); mapped_reads = c()
for (i in 1:length(dirs)) {
    path.to.files <- dirs[i]
    basename <- basename(dirs[i])
    stats.file <- readLines(file.path(path.to.files, paste(basename, "_cutadapt_sorted.stats.txt", sep="")))
    rrna.file <- read.table(file.path(path.to.files, paste(basename, "_cutadapt_sorted.rrnacont.txt", sep="")), check.names = FALSE)
    sample_name <- append(sample_name, basename)
    rrna <- append(rrna, as.integer(rrna.file[1]))
    mapped_reads <- append(mapped_reads, as.integer(strsplit(stats.file[grep(pattern="reads mapped", x = stats.file)[1]], "\t")[[1]][3]))
}

df <- data.frame(sample = sample_name, rRNA = rrna, Total_Mapped = mapped_reads)

png(paste(wd, "/multiqc/rRNAcontamination.png", sep = ""), res=200, unit="in", height=8, width=11)
rrna_plot <- ggplot(df, aes(x=sample_name, y=rRNA/Total_Mapped))+geom_bar(position="stack", stat="identity", color="blue")+
    theme_minimal()+labs(x="Samples", y="Percentage of rRNA reads") + theme(text = element_text(size=6), 
    axis.text.x = element_text(angle=90, hjust=1)) + 
    ggtitle("Contamination Level by rRNA")
rrna_plot + coord_flip()
dev.off()
