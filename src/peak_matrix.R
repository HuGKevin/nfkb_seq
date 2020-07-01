### Script for analyses on downsample experiments.

setwd("/home/kh593/scratch60/nfkb_seq/")

library(dplyr)
library(GenomicRanges)
library(regioneR)
library(data.table)
library(tools)

consensus.beds.dir <- "analysis/overlaps"

index <- fread("/home/kh593/project/nfkb_seq/data/called_libs.tsv",
               col.names = c("donor", "expt", "stim", "lib")) # List of libraries and their metadata
atac_files <- index %>%
    filter(expt == "ATAC") %>%
    filter(!(donor %in% c("TB0611", "TB5728", "TB6578")))
mint_files <- index %>%
    filter(expt == "H3K27ac")

### Do everything for ATAC first
### Count matrix against MCL reference peaks
## Load in a reference peakset
atac_base <- sort(granges(toGRanges(fread(paste0(consensus.beds.dir, "/", atac_files$lib[1], "_count.bed")))))

for(i in 1:nrow(atac_files)){
    print(atac_files$lib[i])
    name <- file_path_sans_ext(atac_files$lib[i])
    
    data <- sort(toGRanges(fread(paste0(consensus.beds.dir, "/", name, "_count.bed"),
                            col.names = c("chr", "start", "end", "width", "strand", "clust_size","count"))))

    mcols(atac_base)[name] <- mcols(data)$count
}

write.table(atac_base, gzfile("analysis/atac_count_matrix.bed.gz"), row.names = FALSE, quote = FALSE)

### Do the same for mintchip
mint_base <- sort(granges(toGRanges(fread(paste0(consensus.beds.dir, "/", mint_files$lib[1], "_count.bed")))))

for(i in 1:nrow(mint_files)){
    print(mint_files$lib[i])
    name <- file_path_sans_ext(mint_files$lib[i])
    
    data <- sort(toGRanges(fread(paste0(consensus.beds.dir, "/", name, "_count.bed"),
                            col.names = c("chr", "start", "end", "width", "strand", "clust_size","count"))))

    mcols(mint_base)[name] <- mcols(data)$count
}

write.table(mint_base, gzfile("analysis/mint_count_matrix.bed"), row.names = FALSE, quote = FALSE)

## ### Correlation between downsampled and consensus peaksets
## ## Load in the data
## metadata <- read.table("/home/kh593/project/downsample/docs/master_index.tsv", header = TRUE, stringsAsFactors = FALSE)[,1:6]
## count.grange <- toGRanges(fread("results/atac_count_matrix.bed"))
## mcols(count.grange)$X21433.threequart <- NULL

## ## Formatting metadata
## metadata$shortcopa <- unlist(lapply(strsplit(metadata$longcopa, split = "_"), "[", 2))

## ## Hierarchical clustering
## count.matrix <- t(as.matrix(mcols(count.grange)[,c(-1,-2)]))
## count.matrix[count.matrix > 0] <- 1

## dist.matrix <- dist(count.matrix, method = "binary")
## clusters <- hclust(dist.matrix)
## plot(clusters)

## cuts <- cutree(clusters, 5)

## copas <- unlist(lapply(strsplit(names(cuts), split = "\\."), "[", 1))
## copas <- substr(copas, 2, nchar(copas))

## downsize <- unlist(lapply(strsplit(names(cuts), split = "\\."), "[", 2))


## cut.frame <- data.frame(cut = cuts,
##                         copa = copas,
##                         size = factor(downsize, levels = c("pooled", "twenfivemil", "tenmil", "fivemil", "onemil", "halfmil",
##                                                            "threequart", "half", "quarter", "tenth")))


## ## Pearson correlation matrix
## cor.matrix <- cor(t(count.matrix))
## cor.matrix[1:10,1:10]




