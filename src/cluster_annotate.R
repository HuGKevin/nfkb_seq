library(DiffBind)
library(vroom)
library(stringr)
library(dplyr)
library(GenomicRanges)
library(ChIPpeakAnno)
library(EnsDb.Hsapiens.v86)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(reactome.db)
library(motifStack)
library(BSgenome.Hsapiens.UCSC.hg38)

setwd("~/scratch60/nfkb_seq/")

peak.dir <- "aligned_reads/"
fig.dir <- "figs/mcl_peak_anno/"

theme_kh <- function(){
  theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

# Run ChipPeakAnno on the filtered atac clusters
## Build gene database
annoData <- toGRanges(EnsDb.Hsapiens.v86, feature = "gene")

## Distribution of mcl peaks around TSS
png(paste0(fig.dir, "atac_tss_enrichment.png"))
binOverFeature(mcl.peaks, annotationData = annoData,
               radius = 5000, nbins = 20, FUN = length, errFun = 0,
               ylab = "count",
               main = "Distribution of mcl peaks around TSS")
dev.off()

## Percentage of mcl peaks around Tx features
aCR <- assignChromosomeRegion(mcl.peaks, nucleotideLevel = FALSE,
                              precedence = c("Promoters", "immediateDownstream", "fiveUTRs", "threeUTRs", "Exons", "Introns"),
                              TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene)
png(paste0(fig.dir, "atac_feature_dist.png"))
barplot(aCR$percentage, las = 3)
dev.off()

## Find nearest promoter in either direction to each feature. 
overlaps.anno <- annotatePeakInBatch(mcl.peaks,
                                     AnnotationData = annoData,
                                     output = "nearestBiDirectionalPromoters",
                                     bindingRegion = c(-2000, 500))
overlaps.anno <- addGeneIDs(overlaps.anno,
                            "org.Hs.eg.db",
                            IDs2Add = "entrez_id")

## 1.4 Proportion of peaks overlapping TSS, or upstream of feature. 
png(paste0(fig.dir, "atac_pie.png"))
pie1(table(overlaps.anno$insideFeature))
dev.off()

## 1.5 enriched GO terms and pathways
over <- getEnrichedGO(overlaps.anno, orgAnn = "org.Hs.eg.db",
                      maxP = 0.05, minGOterm = 10,
                      multiAdjMethod= "BH", condense = TRUE)

path <- getEnrichedPATH(overlaps.anno, "org.Hs.eg.db", "reactome.db", maxP = 0.05)
head(path)
convert <- as.list(reactomePATHID2NAME)
path$PATH <- unlist(convert[as.character(path$path.id)])
significant.paths <- path %>%
    dplyr::select(PATH, pvalue) %>%
    unique %>%
    arrange(pvalue) %>%
    mutate(PATH = factor(PATH))

g1 <- ggplot(data = significant.paths[1:20,], aes(x = -log10(pvalue), y = reorder(PATH, -pvalue))) +
    geom_bar(stat = "identity") +
    theme_kh() +
    labs(title = "20 most significant enriched pathways",
         x = "-log(p)", y = "Pathway")
ggsave(paste0(fig.dir, "atac_top_pathways.png"), g1)

## 1.6 Obtain sequences (Everything below this for the filtered atac clusters doesn't work.)
seq <- getAllPeakSequence(mcl.peaks, upstream = 20, downstream = 20, genome = Hsapiens)

## 1.7 Output a summary of consensus in the peaks
freqs <- oligoFrequency(Hsapiens$chr1, MarkovOrder = 3)
os <- oligoSummary(seq, oligoLength = 6, MarkovOrder = 3,
                   quickMotif = FALSE, freqs = freqs)

pfms <- mapply(function(.ele, id)
    new("pfm", mat = .ele, name = paste("SAMPLE motif", id)),
    os$motifs, 1:length(os$motifs))



# Running ChipPeakAnno on results of DiffBind on ATACseq clusters
library(DiffBind)
load("atac_dba_object.RData")
db <- dba.report(analysis)
fig.dir <- "/figs/dba_figs/"

## PCA on all clusters
png(paste0(fig.dir, "atac_pca.png"))
dba.plotPCA(analysis, DBA_TREATMENT)
dev.off()
## Not much difference between samples on the basis of treatment, actually. Gives support for the idea that accessibility is constitutive.

## PCA on just the differentially accessible clusters
png(paste0(fig.dir, "atac_diff_pca.png"))
dba.plotPCA(analysis, contrast = 1, DBA_TREATMENT)
dev.off()
## Much more pronounced difference in the differentially bound peaks, obviously. Note that PC1, the axis along which treatment seems to vary, explains 41% of the variation.

## MA plot
png(paste0(fig.dir, "atac_maplot.png"))
dba.plotMA(analysis)
dev.off()
## It's apparent the majority of differentially-accessible sites are less accessible in the vehicle condition than in the stimulation condition.
## I think the red dots are the significant differentially-accessible sites and the blue are nonsignificant, after multiple testing correciton.

## Volcano plots
png(paste0(fig.dir, "atac_volcano.png"))
dba.plotVolcano(analysis)
dev.off()
## Same as from the MA plot

## Boxplots
png(paste0(fig.dir, "atac_boxplot.png"))
dba.plotBox(analysis)
dev.off()
## In differentially bound sites, we see there's not much of a difference in mean accessibility, as measured by normalized fragment count, in either condition. The middle columns tell us that the peaks that are more accessible in TNFa conditions are slightly more open in the TNFa condition than the vehicle condition, and vice versa in the rightmost columns.

## Heatmaps
png(paste0(fig.dir, "atac_heatmap.png"))
corvals <- dba.plotHeatmap(analysis)
dev.off()
## I think the lesson here is the same as from the PCA - overall ,not much difference between the treatments.
png(paste0(fig.dir, "atac_diff_heatmap.png"))
corvals <- dba.plotHeatmap(analysis, contrast = 1, correlations = FALSE)
dev.off()
## When we compare the conditions, then we see there are like, striations in the peaks? Interesting. Not sure how to interpret this

## Annotating the differential peaks.
png(paste0(fig.dir, "atac_diff_tss.png"))
binOverFeature(db, annotationData = annoData,
               radius = 5000, nbins = 20, FUN = length, errFun = 0,
               ylab = "count",
               main = "Distribution of differential mcl peaks around TSS")
dev.off()
## Kind of all over the place relative to the TSS.

## Distribution of features under differential clusters
aCR <- assignChromosomeRegion(db, nucleotideLevel = FALSE,
                              precedence = c("Promoters", "immediateDownstream", "fiveUTRs", "threeUTRs", "Exons", "Introns"),
                              TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene)
png(paste0(fig.dir, "atac_diff_barfeat.png"))
barplot(aCR$percentage, las = 3, main = "Feature annotation of differential ATAC-seq peaks")
dev.off()
## Mostly overlapping introns and intergenic regions, but also some promoters.

## Find nearest promoter in either direction within (-5kb, 3kb) or cluster.
db.anno <- annotatePeakInBatch(db,
                               AnnotationData = annoData,
                               output = "nearestBiDirectionalPromoters",
                               bindingRegion = c(-5000, 3000))
db.anno <- addGeneIDs(db.anno,
                      "org.Hs.eg.db",
                      IDs2Add = "entrez_id")

## Plot nearest bidirectional promoter on volcano plot of differential clusters. 
library(ggplot2)
library(ggrepel)

db.df <- as.data.frame(mcols(db.anno))
labels <- as.data.frame(mcols(db.anno)) %>%
                        dplyr::filter(FDR < 0.0001)

g1 <- ggplot(data = as.data.frame(mcols(db)), aes(x = log2(Conc_none/Conc_tnfa), y = -log10(FDR))) +
    geom_point() +
    geom_text_repel(data = labels,
                    aes(label = gene_name),
                    nudge_x = 3.5 - log2(labels$Conc_none/labels$Conc_tnfa),
                    segment.size  = 0.2,
                    segment.color = "grey50",
                    direction     = "y",
                    hjust         = 1) +
    theme_kh() +
    lims(x = c(-.35, .35)) +
    labs(title = "Volcano plot of differential ATACseq peaks, nearest genes for FDR < 10^-4")
ggsave(paste0(fig.dir, "atac_volcano_genes.png"), g1)
