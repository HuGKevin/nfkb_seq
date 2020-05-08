#!/bin/bash

#SBATCH --partition=general
#SBATCH --job-name=pool_reads.%a
#SBATCH --cpus-per-task=16
#SBATCH --mem=20gb
#SBATCH -o /home/kh593/scratch60/nfkb_seq/logs/pool_reads%a.out
#SBATCH -e /home/kh593/scratch60/nfkb_seq/logs/pool_reads%a.err
#SBATCH --array=1-453

module purge
module load SAMtools

index_file="/home/kh593/project/nfkb_seq/data/pooling_array.tsv"

### job information
donor=$(awk -F'\t' -v row=${SLURM_ARRAY_TASK_ID} -v num=1 'FNR == row {print $num}' $index_file)
expt=$(awk -F'\t' -v row=${SLURM_ARRAY_TASK_ID} -v num=2 'FNR == row {print $num}' $index_file)
stim=$(awk -F'\t' -v row=${SLURM_ARRAY_TASK_ID} -v num=3 'FNR == row {print $num}' $index_file)
lib=$(awk -F'\t' -v row=${SLURM_ARRAY_TASK_ID} -v num=4 'FNR == row {print $num}' $index_file)

merged_file="/home/kh593/scratch60/nfkb_seq/pooled_reads/${lib}_pooled.bam"
copas=$(grep "${donor}" /home/kh593/project/nfkb_seq/data/pooling_index.tsv |
	    grep "${expt}" | 
	    grep "${stim}" |
	    cut -f5 |
	    sed "s/$/.final.bam/g" |
	    sed "s/^/\/home\/kh593\/scratch60\/nfkb_seq\/aligned_reads\//g")
echo "${donor} ${expt} ${stim} ${lib}"
echo $copas
echo $merged_file

samtools merge -@ 16 -f ${merged_file} $copas
echo "${lib} completed"
