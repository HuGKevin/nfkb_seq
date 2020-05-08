#!/bin/bash

#SBATCH --partition=general
#SBATCH --job-name=downsample_reads.%a
#SBATCH --cpus-per-task=8
#SBATCH --mem=10gb
#SBATCH -o /home/kh593/scratch60/nfkb_seq/logs/downsample_reads%a.out
#SBATCH -e /home/kh593/scratch60/nfkb_seq/logs/downsample_reads%a.err
#SBATCH --array=1-453

module purge
module load SAMtools

index_file="/home/kh593/project/nfkb_seq/data/pooling_array.tsv"

### job information
donor=$(awk -F'\t' -v row=${SLURM_ARRAY_TASK_ID} -v num=1 'FNR == row {print $num}' $index_file)
expt=$(awk -F'\t' -v row=${SLURM_ARRAY_TASK_ID} -v num=2 'FNR == row {print $num}' $index_file)
stim=$(awk -F'\t' -v row=${SLURM_ARRAY_TASK_ID} -v num=3 'FNR == row {print $num}' $index_file)
lib=$(awk -F'\t' -v row=${SLURM_ARRAY_TASK_ID} -v num=4 'FNR == row {print $num}' $index_file)

target_depth=30000000
fullset="/home/kh593/scratch60/nfkb_seq/pooled_reads/${lib}_pooled.bam"
output="/home/kh593/scratch60/nfkb_seq/downsampled/${lib}.bam"

# compute total reads
total_reads=$(samtools view -@ 8 -c ${fullset})

# compute fraction of reads given an input read depth
frac=$(awk -v down=$target_depth -v full=$total_reads 'BEGIN {frac=down/full;
if (frac > 1) {print 1} else {print frac}}')	  

# samtools view to downsample appropriately
if [ $frac -eq 1 ]
then
    echo "${lib} doesn't exceed downsample threshold"
else
    samtools view -@ 8 -bs $frac $fullset > $output
fi
