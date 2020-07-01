### Script for processing cluster output from MCL. Goal is to produce bed-format clusters.

args <- commandArgs(trailingOnly = TRUE)
expt <- args[1]
chr <- args[2]

setwd(paste0("/home/kh593/scratch60/nfkb_seq/mcl/", expt))

library(IRanges)
library(GenomicRanges)
library(regioneR)
library(data.table)

### Read in data
print(chr)
mcl.output <- readLines(paste0("clusters/", chr, "_clusters.txt"))
og.peaks <- toGRanges(fread(paste0(chr, ".bed")))


peaks <- strsplit(mcl.output, split = "\t")

temp <- unlist(GRangesList(lapply(peaks, function(x){
    subset <- og.peaks[og.peaks$V4 %in% x,]
    result <- GRanges(
        seqnames = unique(seqnames(subset)),
        ranges = IRanges(min(start(subset)):max(end(subset))),
        hits = length(subset))
    return(result)
})))

write.table(temp, paste0(chr, "_CLUSTER.BED"), row.names = FALSE, quote = FALSE, sep = "\t", col.names = FALSE)
