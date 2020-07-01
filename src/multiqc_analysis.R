### We're analyzing the multifastqc data.

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr) 

setwd("/home/kh593/scratch60/nfkb_seq")

theme_kh <- function(){
  theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

data <- fread("qc_stats/summary.filter.tsv")
index <- fread("doc/file_index.csv")

data <- data %>%
    mutate(copa = unlist(map(Library, ~ strsplit(.x, "\\.")[[1]][1]))) %>%
    mutate(copa = unlist(map(copa, ~ substr(.x, 1, nchar(.x)-3)))) %>%
    left_join(index %>%
              dplyr::select(copa, Expt))
    
### Number of pass/fail for each test.
## Sequence Lengths for each library
table(data$SeqLength, data$Expt, data$Trim)

## Average quality at each position of sequence
table(data$PerBaseSeqQual, data$Expt, data$Trim)

## Average quality per tile, for bias at certain positions of the flow cell
table(data$PerTileSeqQual, data$Expt, data$Trim)

## Number of reads with each average sequence quality score
table(data$PerSeqQual, data$Expt, data$Trim)

## Homogeneity of sequences at each position of sequence.
table(data$PerBaseSeq, data$Expt, data$Trim)

## Distirbution of %GC per sequence. Warning or failure if deviation from normality.
table(data$PerSeqGC, data$Expt, data$Trim)

## Distribution of ambiguous calls at each basepair.
table(data$PerBaseN, data$Expt, data$Trim)

## Distribution of sequence lengths, warnings and failures and all.
table(data$SeqLengthDist, data$Expt, data$Trim)

## Distribution of sequence duplication levels
table(data$SeqDupLevel, data$Expt, data$Trim)

## Presence of overrepresented sequences in the first 0.1% of the library
table(data$OverrepSeq, data$Expt, data$Trim)

## Presence of adapter matches in the first 50bp of the first 0.1% of the library drawing from a few databases.
table(data$AdaptCont, data$Expt, data$Trim)

## Histogram of total deduplication percentage, the percentage of the library remaining if duplicated reads are removed.
ggplot(data = dplyr::filter(data, Trim == "pretrim"), aes(x = TotalDedupPerc)) +
    geom_histogram(aes(y = ..density..)) +
    theme_kh() +
    facet_wrap(~ Expt) +
    labs(title = "Sequence Deduplication Distribution",
         x = "% of Reads Remaining After Deduplication")

test <- data %>%
    dplyr::filter(SeqLength %in% c("58", "76"))
