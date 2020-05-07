#!/bin/bash
######################################
### This script is for peak calling our MintChIP experiments. It takes in a paired-end bam file and
### outputs 1) gzipped broadPeaks; 2) gzipped gappedPeaks; 3) signal track of fold-enrichment;
### 4) signal track of signal poisson pvalue. Iterated by peak_call.sh.
### Usage: ./mint_peak_temp.sh 

### Load necessary packages and modules
module load miniconda
source activate atac_env

module load BEDTools
module load MACS2

### Relevant variables and files
## Inputs
input_file="${1}"
COPA_id=$(basename ${input_file} .bam)
pval_thresh="0.01"
genome_sizes="/home/kh593/project/genomes/hg19/hg19.chrom.sizes"

## Intermediate filenames
bpeakfile="${COPA_id}.broadPeak.gz"
gpeakfile="${COPA_id}.gappedPeak.gz"

## Result files
FE_bedgraph="${COPA_id}_FE.bdg"
Pval_bedgraph="${COPA_id}_ppois.bdg"
FE_bigwig="${COPA_id}_FE.bigwig"
Pval_bigwig="${COPA_id}_ppois.bigwig"

### Code that will probably be put into a loop at some point

### Call broad peaks on the data
macs2 callpeak \
      -t $input_file -f BAMPE -n $COPA_id -p $pval_thresh \
      --nomodel --broad --keep-dup all -B --SPMR

### Sort and compress
## Broad
sort -k 8gr,8gr ${COPA_id}_peaks.broadPeak | \
    awk 'BEGIN{OFS = "\t"}{$4 = "Peak_"NR ; print $0}' | \
    gzip -c > $bpeakfile

## Gapped
sort -k 8gr,8gr ${COPA_id}_peaks.gappedPeak | \
    awk 'BEGIN{OFS = "\t"}{$4 = "Peak_"NR ; print $0}' | \
    gzip -c > $gpeakfile

### Convert peaks to signal tracks
## Fold enrichment track
macs2 bdgcmp -t ${COPA_id}_treat_pileup.bdg \
      -c ${COPA_id}_control_lambda.bdg \
      --o-prefix ${COPA_id} -m FE

## Poisson pval track
macs2 bdgcmp -t ${COPA_id}_treat_pileup.bdg \
      -c ${COPA_id}_control_lambda.bdg \
      --o-prefix ${COPA_id} -m ppois

### Sort bedgraph and convert to bigwig
bedSort ${FE_bedgraph} ${FE_bedgraph}
bedSort ${Pval_bedgraph} ${Pval_bedgraph}
bedGraphToBigWig ${FE_bedgraph} $genome_sizes $FE_bigwig
bedGraphToBigWig ${Pval_bedgraph} $genome_sizes $Pval_bigwig

### Remove unnecessary files
rm -f ${COPA_id}_peaks.broadPeak
rm -f ${COPA_id}_peaks.gappedPeak
rm -f ${COPA_id}_peaks.xls
rm -f ${FE_bedgraph}
rm -f ${Pval_bedgraph}
rm -f ${COPA_id}_treat_pileup.bdg
rm -f ${COPA_id}_control_lambda.bdg
