#!/bin/bash

#SBATCH --job-name=bam2bigwig
#SBATCH --mem=10gb --cpus-per-task=1
#SBATCH -o /home/kh593/scratch60/downsample/logs/bam2bigwig.out
#SBATCH -e /home/kh593/scratch60/downsample/logs/bam2bigwig.err

file=$(basename "${1}" .bdg)
bedGraphToBigWig ${1} /home/kh593/project/genomes/hg19/hg19.chrom.sizes "/home/kh593/scratch60/downsample/data/bigwigs/${file}.bw"
