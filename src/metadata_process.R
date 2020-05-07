### Script for making metadata tables for iterating stuff. 
setwd("/home/kh593/project/nfkb_seq/data/")

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

meta <- fread("ALL_METADATA.csv", header = TRUE, check.names = TRUE)

filtered <- meta %>%
  dplyr::select(Donor.External.Collaborator.Identifier.of.Donor, PchAb.Epitope, Reads.1.Filename.URI, BioSam.Cell.Type, Pool.Component.Library) %>%
  mutate(copa = unlist(map(Reads.1.Filename.URI, ~ str_extract(.x, "CoPA_[:digit:]{5}_L[1234]")))) %>%
  dplyr::select(-Reads.1.Filename.URI) %>%
  rename(donor = Donor.External.Collaborator.Identifier.of.Donor, expt = PchAb.Epitope, stim = BioSam.Cell.Type, dnalib = Pool.Component.Library) %>%
  mutate(stim = unlist(map(stim, ~ tail(str_split(.x, " ")[[1]], n = 1)))) %>%
  mutate(dnalib = str_replace(dnalib, " ", "_"))
filtered$stim[filtered$stim == "unstimulated"] <- "none"
filtered$stim[filtered$stim == "stimulated"] <- "tnfa"

write.table(filtered, "/home/kh593/project/nfkb_seq/data/pooling_index.tsv", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

## ### Downsample index
## test <- read.table("/home/kh593/project/downsample2/docs/index.tsv", header = TRUE, sep = "\t")

## test <- test %>%
##   mutate(eighty = 80 * 1000000/reads,
##          seventy = 70 * 1000000/reads,
##          sixty = 60 * 1000000/reads,
##          fifty = 50 * 1000000/reads,
##          forty = 40 * 1000000/reads,
##          thirty = 30 * 1000000/reads,
##          twenty = 20 * 1000000/reads,
##          ten = 10 * 1000000/reads,
##          five = 5 * 1000000/reads,
##          two = 2 * 1000000/reads) %>%
##   dplyr::select(-reads) %>% 
##   gather("milreads", "prop", eighty, seventy, sixty, fifty, forty, thirty, twenty, ten, five, two) %>%
##   filter(prop < 1)

## write.table(test, "/home/kh593/project/downsample2/docs/downsample_index.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

## ### Peak calling index
## test <- read.table("/home/kh593/project/downsample2/docs/index.tsv", header = TRUE, sep = "\t")

## test <- test %>%
##   mutate(eighty = 80 * 1000000/reads,
##          seventy = 70 * 1000000/reads,
##          sixty = 60 * 1000000/reads,
##          fifty = 50 * 1000000/reads,
##          forty = 40 * 1000000/reads,
##          thirty = 30 * 1000000/reads,
##          twenty = 20 * 1000000/reads,
##          ten = 10 * 1000000/reads,
##          five = 5 * 1000000/reads,
##          two = 2 * 1000000/reads) %>%
##   dplyr::rename(full = reads) %>%
##   mutate(full = 0) %>% 
##   gather("milreads", "prop", full, eighty, seventy, sixty, fifty, forty, thirty, twenty, ten, five, two) %>%
##   filter(prop < 1) %>%
##   dplyr::select(-prop)

## write.table(test, "/home/kh593/project/downsample2/docs/peakcall_index.txt", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

