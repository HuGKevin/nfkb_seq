#!/bin/bash
# Usage: sbatch bind_clusters.sh

#SBATCH --job-name=bind_clusters
#SBATCH --partition=scavenge
#SBATCH --mem=25gb --cpus-per-task=1
#SBATCH -e /home/kh593/scratch60/nfkb_seq/logs/bind_clusters%a.err
#SBATCH -o /home/kh593/scratch60/nfkb_seq/logs/bind_clusters%a.out
#SBATCH --array=1-48

module purge

module load R/3.6.1-foss-2018b

index="/home/kh593/project/nfkb_seq/data/mcl_index.tsv"

chr=$(awk -F"\t" -v row=${SLURM_ARRAY_TASK_ID} -v num=1 'FNR == row {print $num}' $index)
expt=$(awk -F"\t" -v row=${SLURM_ARRAY_TASK_ID} -v num=2 'FNR == row {print $num}' $index)

Rscript /home/kh593/project/nfkb_seq/src/bind_clusters.R ${expt} ${chr}

echo "Job completed"
