### Script for analysis of all statistics from FastQC, alignment, and filtering. 
rm(list = ls())

setwd("/home/kh593/scratch60/nfkb_seq")
qc_dir <- "final_qc/"
data_dir <- "/gpfs/ysm/project/cotsapas/kh593/nfkb_seq/data/"
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(purrr)
library(magrittr)

theme_kh <- function(){
  theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

align.stats <- fread(paste0(qc_dir, "bowtie2.alignment.stats.txt"))
picard.stats <- fread(paste0(qc_dir, "picard.alignment.stats.txt"))
qual.per.base <- fread(paste0(qc_dir, "sequence.quality.per.base.txt"))
gc.content <- fread(paste0(qc_dir, "gc.content.level.txt"))
seq.dup <- fread(paste0(qc_dir, "sequence.duplication.level.txt"))
file.index <- fread(paste0(data_dir, "file_index.csv"))
atac.meta <- fread(paste0(data_dir, "atac_lane_metadata.csv"))
mint.meta <- fread(paste0(data_dir, "mint_lane_metadata.csv"))

big.index <- file.index %>%
    dplyr::select(UID, Expt, copa) %>%
    left_join(bind_rows(atac.meta, mint.meta) %>%
              dplyr::select(c(1,2,4,7,8,9,10,14,15,44)))


### Pre-align stats
gc.content
qual.per.base$Base %<>% as.integer
seq.dup

ggplot(data = qual.per.base[qual.per.base$Sample == "CoPA_15577_L1",]) +
    geom_line(aes(x = Base, y = Mean, linetype = Read, color = Phase)) +
    theme_kh()

### Post-align stats
formatted.alignment <- picard.stats %>%
    mutate(blah = unlist(map2(Phase, str_replace_all(Class, " ", "."),
                              ~ paste(.x, .y, sep = ".")))) %>%
    dplyr::select(Sample, Count, blah) %>% 
    pivot_wider(Sample, names_from = blah, values_from = Count) %>%
    left_join(align.stats) %>%
    left_join(big.index, by = c("Sample" = "copa"))

fwrite(formatted.alignment, "results/post_align_stats.csv")

test <- align.stats %>%
    spread(Class, Count) %>%
    dplyr::filter(
