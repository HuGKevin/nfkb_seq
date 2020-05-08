#!/bin/bash

# ### Load necessary packages and modules
# module load miniconda
# source activate atac_env

# module load BEDTools
# module load MACS2

### Relevant variables and files
## Inputs
input_file="${1}"
file_id=$(basename ${input_file} .bam)
pval_thresh="0.01"
genome_sizes="/home/kh593/project/genomes/hg38/hg38_no_alt.chrom.sizes"

## Intermediate filenames
bpeakfile="${file_id}.broadPeaks.gz"
gpeakfile="${file_id}.gappedPeaks.gz"

## Result files
FE_bedgraph="${file_id}_FE.bdg"
Pval_bedgraph="${file_id}_ppois.bdg"
FE_bigwig="${file_id}_FE.bigwig"
Pval_bigwig="${file_id}_ppois.bigwig"

### Code that will probably be put into a loop at some point

### Call broad peaks on the data
macs2 callpeak \
      -t $input_file -f BAMPE -n $file_id -p $pval_thresh \
      --nomodel --broad --keep-dup all -B --SPMR

### Sort and compress
## Broad
sort -k 8gr,8gr ${file_id}_peaks.broadPeak | \
    awk 'BEGIN{OFS = "\t"}{$4 = "Peak_"NR ; print $0}' | \
    gzip -c > $bpeakfile

## Gapped
sort -k 8gr,8gr ${file_id}_peaks.gappedPeak | \
    awk 'BEGIN{OFS = "\t"}{$4 = "Peak_"NR ; print $0}' | \
    gzip -c > $gpeakfile

### Convert peaks to signal tracks
## Fold enrichment track
macs2 bdgcmp -t ${file_id}_treat_pileup.bdg \
      -c ${file_id}_control_lambda.bdg \
      --o-prefix ${file_id} -m FE

## Poisson pval track
macs2 bdgcmp -t ${file_id}_treat_pileup.bdg \
      -c ${file_id}_control_lambda.bdg \
      --o-prefix ${file_id} -m ppois

### Sort bedgraph and convert to bigwig
bedSort ${FE_bedgraph} ${FE_bedgraph}
bedSort ${Pval_bedgraph} ${Pval_bedgraph}
bedGraphToBigWig ${FE_bedgraph} $genome_sizes $FE_bigwig
bedGraphToBigWig ${Pval_bedgraph} $genome_sizes $Pval_bigwig

if [ -f "${FE_bigwig}" ] 
then
    rm -f ${FE_bedgraph}
fi

fi [ -f "${Pval_bigwig}" ]
then
    rm -f ${Pval_bedgraph}
fi

### Remove unnecessary files
rm -f ${file_id}_peaks.broadPeak
rm -f ${file_id}_peaks.gappedPeak
rm -f ${file_id}_peaks.xls
rm -f ${file_id}_treat_pileup.bdg
rm -f ${file_id}_control_lambda.bdg
