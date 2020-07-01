### Script for testing whether sequencing depth of our libraries correlates with biological variables. Goal is to rule out potential biases.
setwd("/home/kh593/project/nfkb_seq/")

library(dplyr)
library(data.table)

atac.depth <- fread("data/atac_depths.csv")
mint.depth <- fread("data/mint_depths.csv")

bio.data <- fread("/home/kh593/project/flow_norm/data/raw/donor_demographics.csv") %>%
    dplyr::select(-subjID) %>%
    mutate(donor = paste0(substr(donorID, 1, 2),
                          substr(donorID, 7, 10))) %>%
    dplyr::select(-donorID)


metadata <- fread("data/ALL_METADATA.csv") %>%
    dplyr::select(donor = `Donor:External Collaborator Identifier of Donor`,
                  expt = `PchAb:Epitope`,
                  stim = `BioSam:Cell Type`,
                  lib = `Pool Component:Library`) %>%
    mutate(stim = unlist(map(stim, ~ lapply(str_split(.x, " "), "[", 5)))) %>%
    mutate(lib = str_replace(lib, " ", "_")) %>%
    unique
metadata[metadata$donor == "TB0611", "donor"] <- "TB6110"
metadata[metadata$donor == "TB7586", "donor"] <- "TB7556"
metadata[metadata$donor == "TB0114", "donor"] <- "TB1114"
metadata[metadata$donor == "TB4509", "donor"] <- "TB5409"

t1 <- left_join(atac.depth, metadata)
t2 <- left_join(t1, bio.data)

test1 <- aov(depth ~ rs1800693, data = t2)
### P = 0.0135 for rs1800693, after correcting mislabeled donors, P = 0.00492

ggplot(data = t2, aes(x = rs1800693, y = depth)) +
    geom_boxplot() +
    scale_y_continuous(labels = comma)

k1 <- kruskal.test(depth ~ rs1800693, data = t2)

## What else should I want to test?
                                        # Seq depth ~ peaks called in library
                                        # baseline nfkb
                                        # delta_nfkb
