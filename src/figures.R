## Workspace for making any and all figures from sequencing data.

setwd("/home/kh593/scratch60/nfkb_seq/")

## Relevant libraries
library(ggplot2)
library(ggcorrplot)
library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(magrittr)
library(regioneR)
library(scales)
library(stringr)
library(logisticPCA)

theme_kh <- function(){
  theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

fig.dir <- "/home/kh593/scratch60/nfkb_seq/figs/"

# Distribution of clusterwidthsd
atac_clusters <- fread("mcl/final_atac.bed")
colnames(atac_clusters) <- c("chr", "start", "end", "width", "strand", "peaks")

g1 <- ggplot(data = atac_clusters, aes(x = width)) +
    geom_histogram() +
    geom_boxplot(aes(y = 250000)) + 
    theme_kh() +
    scale_y_continuous(labels = comma) +
    labs(title = "Clusterwidths of ATAC-seq clusters, prefiltering",
         x = "Width (bp)", y = "Count")
ggsave(paste0(fig.dir, "atac_clusterwidths.png"), g1)

# How many clusters are present in x samples? 
atac_active <- toGRanges(fread("analysis/atac_count_matrix.bed.gz"))
presence <- data.frame(granges(atac_active),
                       class = apply(as.data.frame(mcols(atac_active)[,3:164]) > 0, 1, sum))

ggplot(data = presence, aes(x = class)) +
    geom_histogram(binwidth = 1) +
    theme_kh() +
    scale_y_log10(labels = comma) +
    labs(title = "Distribution of cluster classes across all libraries")

# What if we do just look at a single library? How many of clusters of each class are present in the library?
compare <- data.frame(presence, real = atac_active$DNA_Lib_10973 > 0)
plotter <- as.data.frame(table(compare$class, compare$real))
colnames(plotter) <- c("class", "present", "freq")
ggplot(data = plotter, aes(x = class, y = freq, fill = present)) +
    geom_bar(stat = "identity", position = "stack") +
    theme_kh() +
    labs(title = "Distribution of cluster classes in DNA_Lib_10973")

# What if we do the above, but plot the proportion of clusters that are present, rather than
# absolute number of clusters? 
p2 <- plotter %>%
    spread(present, freq) %>%
    dplyr::rename("pres" = "TRUE", "notpres" = "FALSE") %>%
    mutate(prop = pres / (pres + notpres))
ggplot(data = p2, aes(x = class, y = prop)) +
    geom_bar(stat = "identity") +
    theme_kh() +
    labs(title = "Proportion of clusters of each class present in DNA_Lib_10973")

# Do PCA on cluster binding profiles
index <- fread("/home/kh593/project/nfkb_seq/data/atac_libs.tsv",
               col.names = c("donor", "expt", "stim", "lib")) # List of libraries and their metadata

count.matrix <- as_tibble(atac_active) %>%
    dplyr::select(-seqnames, -start, -end, -width, -strand, -width.1, -strand.1) %>% 
    t
count.matrix <- as.data.frame(count.matrix != 0 )
count.matrix <- rownames_to_column(count.matrix, "lib")
count.matrix <- bind_cols(count.matrix, index)
count.matrix <- count.matrix %>% filter(str_detect(donor, "TB"))

pc1 <- prcomp(count.matrix[,2:2141433])
plotter <- bind_cols(count.matrix[,c("donor", "stim", "lib")], as.data.frame(pc1$x))

vars <- pc1$sdev^2 ## For getting proportion of variance explained
scree <- cumsum(vars) / sum(vars)
propvar <- vars / sum(vars)

p1.12 <- ggplot(data = plotter) +
    geom_point(aes(x = PC1, y = PC2, color = stim)) +
    labs(x = paste0("PC1: ", round(100 * propvar[1], 2), "%"),
         y = paste0("PC2: ", round(100 * propvar[2], 2), "%"),
         title = "PCA of cluster replication in downsamples") + 
    theme_kh() # Nothing looks reasonable here.
    
p1.23 <- ggplot(data = plotter) +
    geom_point(aes(x = PC2, y = PC3, color = stim)) +
    labs(x = paste0("PC2: ", round(100 * propvar[2], 2), "%"),
         y = paste0("PC3: ", round(100 * propvar[3], 2), "%"),
         title = "PCA of cluster replication in downsamples") + 
    theme_kh() # Nothing looks reasonable here.

p1.34 <- ggplot(data = plotter) +
    geom_point(aes(x = PC3, y = PC4, color = stim)) +
    labs(x = paste0("PC3: ", round(100 * propvar[3], 2), "%"),
         y = paste0("PC4: ", round(100 * propvar[4], 2), "%"),
         title = "PCA of cluster replication in downsamples") + 
    theme_kh() # Nothing looks reasonable here.

p1.45 <- ggplot(data = plotter) +
    geom_point(aes(x = PC4, y = PC5, color = stim)) +
    labs(x = paste0("PC4: ", round(100 * propvar[4], 2), "%"),
         y = paste0("PC5: ", round(100 * propvar[5], 2), "%"),
         title = "PCA of cluster replication in downsamples") + 
    theme_kh() # Nothing looks reasonable here.

scree_plotter <- data.frame(cumpropvar = scree, propvar = propvar, pc = 1:length(scree))
sp1 <- ggplot(data = scree_plotter, aes(x = pc, y = cumpropvar)) +
    geom_point() +
    geom_line() +
    theme_kh() +
    labs(title = "Screeplot of PCA on ATAC-clusters") 

ggsave("figures/pc1_12.png", p1.12)
ggsave("figures/pc1_23.png", p1.23)
ggsave("figures/pc1_34.png", p1.34)
ggsave("figures/pc1_45.png", p1.45)
ggsave("figures/pc1_scree.png", sp1)

rm(pc1)
### MY questions are
### 1) Without filtering, the first PC seems like it's just those techdev samples. What happens if i look at lower PCs?
###### Going down to 3x4, nothing interesting really appears. I think other techdev samples are still killing it.
###### Each PC comparison just results in most of the libraries in the middle, and a few on the sides. So that doesn't seem like it's right. 
###### I should look at the pivot table at dna_lib by date run.

### 2) What if I do that while removing apparent techdev samples?
index <- fread("/home/kh593/project/nfkb_seq/data/peakcall_array.tsv",
               col.names = c("donor", "expt", "stim", "lib")) # List of libraries and their metadata

count.matrix <- atac_active %>%
    dplyr::select(-seqnames, -start, -end, -width, -strand) %>% 
    t
count.matrix <- as.data.frame(count.matrix != 0 )
count.matrix <- rownames_to_column(count.matrix, "lib")
count.matrix <- bind_cols(count.matrix, index[index$expt == "ATAC",])
count.matrix <- count.matrix %>% filter(str_detect(donor, "TB"))

pc2 <- prcomp(count.matrix[,2:2141433])
plotter <- bind_cols(count.matrix[,c("donor", "stim", "lib")], as.data.frame(pc2$x))

vars <- pc2$sdev^2 ## For getting proportion of variance explained
scree <- cumsum(vars) / sum(vars)
propvar <- vars / sum(vars)

p2.12 <- ggplot(data = plotter) +
    geom_point(aes(x = PC1, y = PC2, color = stim)) +
    labs(x = paste0("PC1: ", round(100 * propvar[1], 2), "%"),
         y = paste0("PC2: ", round(100 * propvar[2], 2), "%"),
         title = "PCA of cluster replication in downsamples, techdevs removed") + 
    theme_kh() # still, nothing makes sense
p2.23 <- ggplot(data = plotter) +
    geom_point(aes(x = PC2, y = PC3, color = stim)) +
    labs(x = paste0("PC2: ", round(100 * propvar[2], 2), "%"),
         y = paste0("PC3: ", round(100 * propvar[3], 2), "%"),
         title = "PCA of cluster replication in downsamples, techdevs removed") + 
    theme_kh() # still, nothing makes sense
p2.34 <- ggplot(data = plotter) +
    geom_point(aes(x = PC3, y = PC4, color = stim)) +
    labs(x = paste0("PC3: ", round(100 * propvar[3], 2), "%"),
         y = paste0("PC4: ", round(100 * propvar[4], 2), "%"),
         title = "PCA of cluster replication in downsamples, techdevs removed") + 
    theme_kh() # still, nothing makes sense

scree_plotter <- data.frame(cumpropvar = scree, propvar = propvar, pc = 1:length(scree))
sp2 <- ggplot(data = scree_plotter, aes(x = pc, y = cumpropvar)) +
    geom_point() +
    geom_line() +
    theme_kh() +
    lab(title = "Screeplot of PCA on ATAC-clusters, techdevs removed") 

ggsave("figures/pc2_12.png", p2.12)
ggsave("figures/pc2_23.png", p2.23)
ggsave("figures/pc2_34.png", p2.34)
ggsave("figures/pc2_scree.png", sp2)

rm(pc2)
### We also just don't see any meaningful spread. Everything is clumped together. The first PCs don't explain very much.
### Suspicious.
### Maybe the issue is I shouldn't be using PCA?
### logisticPCA, look into that. 

### 2b) What if we filter everything that is before nov-1-19? I think after that is when everything is DEFINITELY not techdev. 
index <- fread("/home/kh593/project/nfkb_seq/data/peakcall_array.tsv",
               col.names = c("donor", "expt", "stim", "lib")) # List of libraries and their metadata

count.matrix <- atac_active %>%
    dplyr::select(-seqnames, -start, -end, -width, -strand) %>% 
    t
count.matrix <- as.data.frame(count.matrix != 0 )
count.matrix <- rownames_to_column(count.matrix, "lib")
count.matrix <- bind_cols(count.matrix, index[index$expt == "ATAC",])
count.matrix <- count.matrix %>%
    filter(str_detect(donor, "TB")) %>%
    filter(lib %in% test[date > as.Date("2019-11-1"), lib])

pc3 <- prcomp(count.matrix[,2:2141433])
plotter <- bind_cols(count.matrix[,c("donor", "stim", "lib")], as.data.frame(pc3$x))

vars <- pc3$sdev^2 ## For getting proportion of variance explained
scree <- cumsum(vars) / sum(vars)
propvar <- vars / sum(vars)

p3.12 <- ggplot(data = plotter) +
    geom_point(aes(x = PC1, y = PC2, color = stim)) +
    labs(x = paste0("PC1: ", round(100 * propvar[1], 2), "%"),
         y = paste0("PC2: ", round(100 * propvar[2], 2), "%"),
         title = "PCA of cluster replication in downsamples") + 
    theme_kh() # things still look weird. 

### 3) What if I just remove all samples that have more than one dnalib per condition?
###  This hasn't worked yet. It's a mess. 
metadata <- fread("/home/kh593/project/nfkb_seq/data/ALL_METADATA.csv", check.names = TRUE)
test <- unique(metadata[,c("Created.At", "PchAb.Epitope", "Donor.External.Collaborator.Identifier.of.Donor",
                           "Pool.Component.Library", "BioSam.Cell.Type")])
colnames(test) <- c("date", "expt", "donor", "lib", "stim")
test$stim[grepl("unstimulated", test$stim)] <- "none"
test$stim[grepl("stimulated", test$stim)] <- "tnfa"
test$lib <- gsub(" ", "_", test$lib)
test$date <- as.Date(test$date, "%m/%d/%y")

test <- test %>%
    filter(expt == "ATAC") %>% 
    group_by(donor, stim) %>% nest() %>%
    mutate(unique_libs = unlist(map(data, ~ length(unique(.x[,"lib"])))))

### 4) What if I filter MCL peaks with less than 4 samples active?
index <- fread("/home/kh593/project/nfkb_seq/data/peakcall_array.tsv",
               col.names = c("donor", "expt", "stim", "lib")) # List of libraries and their metadata

## Load some additional metadata for date filtering. 
metadata <- fread("/home/kh593/project/nfkb_seq/data/ALL_METADATA.csv", check.names = TRUE)
test <- unique(metadata[,c("Created.At", "PchAb.Epitope", "Donor.External.Collaborator.Identifier.of.Donor",
                           "Pool.Component.Library", "BioSam.Cell.Type")])
colnames(test) <- c("date", "expt", "donor", "lib", "stim")
test$stim[grepl("unstimulated", test$stim)] <- "none"
test$stim[grepl("stimulated", test$stim)] <- "tnfa"
test$lib <- gsub(" ", "_", test$lib)
test$date <- as.Date(test$date, "%m/%d/%y")

count.matrix <- atac_active %>%
    dplyr::select(-seqnames, -start, -end, -width, -strand) %>% 
    t
count.matrix <- as.data.frame(count.matrix != 0 )
count.matrix <- rownames_to_column(count.matrix, "lib")
count.matrix <- bind_cols(count.matrix, index[index$expt == "ATAC",])
count.matrix <- count.matrix %>%
    filter(str_detect(donor, "TB")) %>%
    filter(lib %in% test[date > as.Date("2019-11-1"), lib])

t <- apply(count.matrix[,2:2141433], 2, sum)
count.matrix <- bind_cols(count.matrix[,c("donor", "lib", "stim")], count.matrix[,which(t > 4) + 1])

pc4 <- prcomp(count.matrix[,c(-1, -2, -3)])
plotter <- bind_cols(count.matrix[,c("donor", "stim", "lib")], as.data.frame(pc4$x))

vars <- pc4$sdev^2 ## For getting proportion of variance explained
scree <- cumsum(vars) / sum(vars)
propvar <- vars / sum(vars)

### dataframe for adding stim arrows
plotter2 <- plotter %>%
    gather("PC", "value", -donor, -stim, -lib) %>%
    mutate(pcstim = paste(PC, stim, sep = ".")) %>%
    dplyr::select(donor, value, pcstim) %>%
    spread(pcstim, value) 

p2.12 <- ggplot(data = plotter) +
    geom_point(aes(x = PC1, y = PC2, color = stim)) +
    labs(x = paste0("PC1: ", round(100 * propvar[1], 2), "%"),
         y = paste0("PC2: ", round(100 * propvar[2], 2), "%"),
         title = "PCA of cluster replication in downsamples, clusters filtered") + 
    theme_kh() +  # PC1 separates well by library size
    geom_segment(data = plotter2,
                 aes(x = PC1.none, y = PC2.none,
                     xend = PC1.tnfa, yend = PC2.tnfa),
                 arrow = arrow(length = unit(0.03, "npc")))

p2.23 <- ggplot(data = plotter) +
    geom_point(aes(x = PC2, y = PC3, color = stim)) +
    labs(x = paste0("PC2: ", round(100 * propvar[2], 2), "%"),
         y = paste0("PC3: ", round(100 * propvar[3], 2), "%"),
         title = "PCA of cluster replication in downsamples, clusters filtered") + 
    theme_kh() +  # PC1 separates well by library size
    geom_segment(data = plotter2,
                 aes(x = PC2.none, y = PC3.none,
                     xend = PC2.tnfa, yend = PC3.tnfa),
                 arrow = arrow(length = unit(0.03, "npc")))

scree_plotter <- data.frame(cumpropvar = scree, propvar = propvar, pc = 1:length(scree))
sp2 <- ggplot(data = scree_plotter, aes(x = pc, y = cumpropvar)) +
    geom_point() +
    geom_line() +
    theme_kh()

ggsave("figures/pc2_12.png", p2.12)
ggsave("figures/pc2_23.png", p2.23)
ggsave("figures/pc2_scree.png", sp2)

# Setup for some more plots below
bigdata <- fread("data/ALL_METADATA.csv", check.names = TRUE)
bigdata <- unique(bigdata[,c("Created.At", "PchAb.Epitope", "Donor.External.Collaborator.Identifier.of.Donor",
                           "Pool.Component.Library", "BioSam.Cell.Type")])
colnames(bigdata) <- c("date", "expt", "donor", "lib", "stim")
bigdata$stim[grepl("unstimulated", bigdata$stim)] <- "none"
bigdata$stim[grepl("stimulated", bigdata$stim)] <- "tnfa"
bigdata$lib <- gsub(" ", "_", bigdata$lib)
bigdata$date <- as.Date(bigdata$date, "%m/%d/%y")

atac.peaks <- gsub(".narrowPeaks", "", list.files("/home/kh593/scratch60/nfkb_seq/results/peak_call/beds/atac/", pattern = "*.narrowPeaks"))
mint.peaks <- gsub(".broadPeaks", "", list.files("/home/kh593/scratch60/nfkb_seq/results/peak_call/beds/mint/", pattern = "*.broadPeaks"))

atac.survive <- bigdata %>%
    filter(expt == "ATAC",
           lib %in% atac.peaks,
           str_detect(donor, "TB")) %>%
    dplyr::select(-date) %>%
    unique
nrow(atac.survive)
length(unique(atac.survive$donor))

mint.survive <- bigdata %>%
    filter(expt == "H3K27ac",
           lib %in% mint.peaks,
           str_detect(donor, "TB")) %>%
    dplyr::select(-date) %>%
    unique
nrow(mint.survive)
length(unique(mint.survive$donor))

# Number of peaks per library
test <- fread("/home/kh593/scratch60/nfkb_seq/results/peak_call/beds/atac/libsize_dist.txt", col.names = c("peaks", "lib"))
test$lib <- gsub(".narrowPeaks", "", test$lib)
atac.reads <- left_join(atac.survive, test)

test <- fread("/home/kh593/scratch60/nfkb_seq/results/peak_call/beds/mint/libsize_dist.txt", col.names = c("peaks", "lib"))
test$lib <- gsub(".broadPeaks", "", test$lib)
mint.reads <- left_join(mint.survive, test)

p1 <- ggplot(data = atac.reads, aes(x = peaks)) +
    geom_boxplot() +
    geom_dotplot(stackdir = "center", dotsize = 0.2) + 
    theme_kh() +
    labs(title = "Number of ATAC peaks per library",
         x = "Number of Peaks", y = "") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.x = element_line())
ggsave(paste0(fig.dir, "atac_peaknum.png"), p1, width = 8, height = 2, units = "in")

p2 <- ggplot(data = mint.reads, aes(x = peaks)) +
    geom_boxplot() +
    geom_dotplot(stackdir = "center", dotsize = 0.2) + 
    theme_kh() +
    labs(title = "Number of Mint peaks per library",
         x = "Number of Peaks", y = "") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.x = element_line()) +
    scale_x_continuous(labels = comma)
ggsave(paste0(fig.dir, "mint_peaknum.png"), p2, width = 8, height = 2, units = "in")

# Plots
atac.cov <- fread("results/analysis/atac_genome_coverage.tsv")
mint.cov <- fread("results/analysis/mint_genome_coverage.tsv")

atac.reads <- left_join(atac.reads, atac.cov)
mint.reads <- left_join(mint.reads, mint.cov)

c1 <- ggplot(data = atac.reads, aes(x = cov)) +
    geom_boxplot() +
    geom_dotplot(stackdir = "center", dotsize = 0.2, binwidth = 0.001) +
    theme_kh() +
    labs(title = "Genome coverage proportion of ATAC peaks per library",
         x = "Proportion of genome covered") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.x = element_line()) +
    xlim(c(0,0.03))
ggsave(paste0(fig.dir, "atac_genome_cov.png"), c1, width = 8, height = 2, units = "in")

c2 <- ggplot(data = mint.reads, aes(x = cov)) +
    geom_boxplot() +
    geom_dotplot(stackdir = "center", dotsize = 0.2) +
    theme_kh() +
    labs(title = "Genome coverage proportion of Mint peaks per library",
         x = "Proportion of genome covered") +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid.major.x = element_line()) +
    xlim(c(0, 0.4))
ggsave(paste0(fig.dir, "mint_genome_cov.png"), c2, width = 8, height = 2, units = "in")

# Cluster data
atac.clust.mat <- fread("results/analysis/atac_count_matrix.bed.gz")

test <- atac.clust.mat %>%
    dplyr::select(-seqnames, -start, -end, -width, -strand)
test2 <- apply(test, 2, function(x) sum(x > 0))
atac.clust.per.lib <- data.frame(clusts = test2) %>%
    rownames_to_column %>%
    dplyr::rename("lib" = "rowname")
atac.reads <- left_join(atac.reads, atac.clust.per.lib)

ac1 <- ggplot(data = atac.reads, aes(x = clusts)) +
    geom_histogram() +
    theme_kh() +
    labs(title = "Distribution of active clusters per ATAC library",
         x = "Clusters", y = "Count (of libraries)")
ggsave(paste0(fig.dir, "atac_clust_per_lib.png"), ac1)

atac.libs.per.clust <- apply(select(test,atac.reads$lib), 1, function(x) sum(x > 0))
lc1 <- ggplot(data = as.data.frame(atac.libs.per.clust), aes(x = atac.libs.per.clust)) +
    geom_histogram() +
    theme_kh() +
    labs(title = "Distribution of number of libraries in which each cluster is active",
         x = "Libraries", y = "Count (of clusters)")
ggsave(paste0(fig.dir, "atac_libs_per_clust.png"), lc1)

# Cluster data
mint.clust.mat <- fread("results/analysis/mint_count_matrix.bed.gz")

test <- mint.clust.mat %>%
    dplyr::select(-seqnames, -start, -end, -width, -strand)
test2 <- apply(test, 2, function(x) sum(x > 0))
mint.clust.per.lib <- data.frame(clusts = test2) %>%
    rownames_to_column %>%
    dplyr::rename("lib" = "rowname")
mint.reads <- left_join(mint.reads, mint.clust.per.lib)

ac2 <- ggplot(data = mint.reads, aes(x = clusts)) +
    geom_histogram() +
    theme_kh() +
    labs(title = "Distribution of active clusters per MINT library",
         x = "Clusters", y = "Count (of libraries)") +
    scale_x_continuous(labels = comma)
ggsave(paste0(fig.dir, "mint_clust_per_lib.png"), ac2)

mint.libs.per.clust <- apply(select(test,mint.reads$lib), 1, function(x) sum(x > 0))
lc2 <- ggplot(data = as.data.frame(mint.libs.per.clust), aes(x = mint.libs.per.clust)) +
    geom_histogram() +
    theme_kh() +
    labs(title = "Distribution of number of libraries in which each cluster is active",
         x = "Libraries", y = "Count (of clusters)")
ggsave(paste0(fig.dir, "mint_libs_per_clust.png"), lc2)



#### Nothing below this ever worked out because it was too computationally expensive. Maybe ya'll can make it work? 
#### TRYING out logisticPCA
atac_prof <- fread("analysis/atac_count_matrix.bed.gz")
index <- fread("/home/kh593/project/nfkb_seq/data/peakcall_array.tsv",
               col.names = c("donor", "expt", "stim", "lib")) # List of libraries and their metadata
genos <- fread("/home/kh593/project/flow_norm/data/NFKB_cells_final.txt")
geno <- genos %>%
    dplyr::select(donor, rs179195, rs1800693) %>%
    unique %>%
    mutate(donorf = sprintf("TB%04d", donor)) %>%
    dplyr::select(donorf, rs179195, rs1800693) %>%
    dplyr::rename(donor = donorf)

index <- left_join(index, geno)

count.matrix <- atac_prof %>%
    dplyr::select(-seqnames, -start, -end, -width, -strand) %>% 
    t
count.matrix <- as.data.frame(count.matrix != 0 )
count.matrix <- rownames_to_column(count.matrix, "lib")
count.matrix <- bind_cols(count.matrix, index[index$expt == "ATAC",])
final.mat <- count.matrix %>%
    filter(str_detect(donor, "TB")) %>%
    filter(!(donor %in% c("TB0611", "TB5728", "TB6578")))

test <- final.mat[,which(colSums(final.mat[,2:2141433]) > 2) + 1]
test2 <- bind_cols(final.mat[,c("lib", "stim", "donor", "expt")], test)
table(apply(final.mat[,2:2141433], 2, function(x) sum(x)))

m1 <- logisticSVD(test2[,5:664680], k = 2)

save(m1, file = "logsvdobj")

plot(m1, type = "scores") + geom_point(aes(color = final.mat$rs1800693)) + ggtitle("Exponential Family PCA on ATAC Cluster Presence")
### Trying out the logisticpca method doesn't work. Simply too much data, ends up producing a vector of size 4TB.

cdist <- dist(test, method = "binary")
hc <- hclust(cdist, method = "complete")
dc <- as.dendrogram(hc)

plot(hc, cex = 0.6, hang = -1)
