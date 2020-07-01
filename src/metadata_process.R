### Script for making metadata tables for iterating stuff. 
setwd("/home/kh593/project/nfkb_seq/data/")

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

meta <- fread("ALL_METADATA.csv", header = TRUE, check.names = TRUE)

base <- meta %>%
    dplyr::select(Donor.External.Collaborator.Identifier.of.Donor, PchAb.Epitope, Reads.1.Filename.URI, BioSam.Cell.Type, Pool.Component.Library) %>%
    mutate(copa = unlist(map(Reads.1.Filename.URI, ~ str_extract(.x, "CoPA_[:digit:]{5}_L[1234]")))) %>%
    dplyr::select(-Reads.1.Filename.URI) %>%
    rename(donor = Donor.External.Collaborator.Identifier.of.Donor,
           expt = PchAb.Epitope, stim = BioSam.Cell.Type, dnalib = Pool.Component.Library) %>%
    mutate(stim = unlist(map(stim, ~ lapply(str_split(.x, " "), "[", 5)))) %>%

    mutate(dnalib = str_replace(dnalib, " ", "_"))
base$stim[base$stim == "unstimulated"] <- "none"
base$stim[base$stim == "TNFa"] <- "tnfa"
base$stim[base$stim == "PMA"] <- "pma"

atac.copa <- base %>%
    filter(expt == "ATAC") %>%
    filter(str_detect(donor, "TB"))

write.table(atac.copa, "/home/kh593/project/nfkb_seq/data/atac_copa.tsv", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

mint.copa <- base %>%
    filter(expt == "H3K27ac") %>%
    filter(str_detect(donor, "TB"))

write.table(mint.copa, "/home/kh593/project/nfkb_seq/data/mint_copa.tsv", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

########## ARRAY FOR POOLING AND DOWNSAMPLING ##################
atac.libs <- unique(dplyr::select(atac.copa, -copa))
mint.libs <- unique(dplyr::select(mint.copa, -copa))

write.table(atac.libs, "/home/kh593/project/nfkb_seq/data/atac_libs.tsv", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
write.table(mint.libs, "/home/kh593/project/nfkb_seq/data/mint_libs.tsv", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

########## ALIGNING, ALIGNMENT QC, PEAK CALLING  ###############
cluster_input <- owjfowjief

atac.downed.libs
mint.downed.libs
