setwd("/home/kh593/project/nfkb_seq/")

library(data.table)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(GenomicRanges)
library(regioneR)
library(scales) 

theme_kh <- function(){
  theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

# Compute proportion of genome covered by each library
## ATAC
atac.peaks <- gsub(".narrowPeaks", "", list.files("/home/kh593/scratch60/nfkb_seq/results/peak_call/beds/atac/", pattern = "*.narrowPeaks"))
genome <- fread("/home/kh593/project/genomes/hg38/hg38_no_alt.chrom.sizes")
genomesize <- sum(genome[1:24,"V2"])
pct <- c()

for(i in 1:length(atac.peaks)){ 
    temp <- toGRanges(fread(paste0("/home/kh593/scratch60/nfkb_seq/results/peak_call/beds/atac/", atac.peaks[i], ".narrowPeaks")))
    coverage <- sum(end(temp) - start(temp))
    pct[i] <- coverage/genomesize
    print(atac.peaks[i])
}

atac.cov <- data.frame(lib = atac.peaks, cov = pct)

write.table(atac.cov, "/home/kh593/scratch60/nfkb_seq/results/analysis/atac_genome_coverage.tsv", sep = "\t", row.names = FALSE)

## Mint
mint.peaks <- gsub(".broadPeaks", "", list.files("/home/kh593/scratch60/nfkb_seq/results/peak_call/beds/mint/", pattern = "*.broadPeaks"))
genome <- fread("/home/kh593/project/genomes/hg38/hg38_no_alt.chrom.sizes")
genomesize <- sum(genome[1:24,"V2"])
pct <- c()

for(i in 1:length(mint.peaks)){ 
    temp <- toGRanges(fread(paste0("/home/kh593/scratch60/nfkb_seq/results/peak_call/beds/mint/", mint.peaks[i], ".broadPeaks")))
    coverage <- sum(end(temp) - start(temp))
    pct[i] <- coverage/genomesize
    print(mint.peaks[i])
}

mint.cov <- data.frame(lib = mint.peaks, cov = pct)

write.table(mint.cov, "/home/kh593/scratch60/nfkb_seq/results/analysis/mint_genome_coverage.tsv", sep = "\t", row.names = FALSE)
