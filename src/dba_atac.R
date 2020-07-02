library(DiffBind)
library(vroom)
library(stringr)
library(dplyr)
library(GenomicRanges)

setwd("~/scratch60/nfkb_seq/")

peak.dir <- "aligned_reads/"
fig.dir <- "figs/mcl_peak_anno/"

theme_kh <- function(){
  theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

# Load metadata and coerce into a data frame for use with DiffBind.
libs <- unlist(lapply(str_split(list.files(peak.dir, pattern = "^DNA_Lib_(.*).down.bam$"), "\\."), "[", 1))
mcl.peaks <- as(vroom("mcl/filtered_atac_clusters.bed", col_names = c("chr", "start", "end", "class")), "GRanges")

metadata <- vroom("~/project/nfkb_seq/data/atac_libs.tsv") %>%
    dplyr::filter(dnalib %in% libs) %>%
    dplyr::rename(Treatment = stim, Factor = expt) %>% 
    mutate(bamReads = paste0(peak.dir, dnalib, ".down.bam")) %>%
    mutate(Peaks = paste0("peaks/atac/", dnalib, ".narrowPeaks")) %>%
    mutate(PeakCaller = "narrow") %>%
    dplyr::rename(SampleID = dnalib) %>%
    dplyr::filter(!(donor %in% c("TB0611", "TB5728", "TB6578")))

# DiffBind analysis
## Create DiffBind object
object <- dba(sampleSheet = metadata)

## Count reads within supplied peakset
object.counted <- dba.count(object, peaks = mcl.peaks)

## Establish contrast matrix - contrast along treatment conditions.
object <- dba.contrast(object.counted, categories = DBA_TREATMENT)

## Perform DeSeq2 analysis.
analysis <- dba.analyze(object)

## Extract differentially accessible peaks.
db <- dba.report(analysis)

## Saving objects because the above analysis required a long time and 100GB+ of memory.
save(object.counted, analysis, db, "atac_dba_object.RData")
write.table(db, "atac_dba_results.tsv",
            quotes = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

