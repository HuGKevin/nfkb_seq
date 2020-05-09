#!/bin/bash

#SBATCH --partition=general
#SBATCH --job-name=mcl_%a
#SBATCH --cpus-per-task=10
#SBATCH --mem=24gb
#SBATCH -o /home/kh593/scratch60/nfkb_seq/logs/mcl_atac%a.out
#SBATCH -e /home/kh593/scratch60/nfkb_seq/logs/mcl_atac%a.err
#SBATCH --array=1-24

module purge
module load MCL
module load BEDTools

index_file="/home/kh593/project/nfkb_seq/data/peakcall_array.tsv"

if [ $SLURM_ARRAY_TASK_ID -eq 23 ]
then
    chr=X
elif [ $SLURM_ARRAY_TASK_ID -eq 24 ]
then
    chr=Y
else
    chr=$SLURM_ARRAY_TASK_ID
fi

mcl_dir="/home/kh593/scratch60/nfkb_seq/results/MCL"
file="${mcl_dir}/atac/chr${chr}.bed"
trimmed="${mcl_dir}/atac/chr${chr}_trimmed.bed"
intersection="${mcl_dir}/atac/chr${chr}_int.bed"
abc="${mcl_dir}/atac/chr${chr}_input.abc"
mcl_output="${mcl_dir}/atac/out.chr${chr}_input.abc.I20"
final="${mcl_dir}/atac/chr${chr}_clusters.txt"

### Cut out excess columns
cut -f1,2,3,4 $file > ${trimmed}

### Intersect against self and report overlaps and size of overlap
intersectBed -a ${trimmed} -b ${trimmed} -wo -sorted > ${intersection}

### Remove unnecessary columns
cut -f4,8,9 ${intersection} > ${abc}

### Run MCL (for future note, can use -te for multithreading)
mcl ${abc} --abc -o ${final} -te 10 

### Rename files
# mv ${mcl_output} ${final} 

# ### Remove intermediate files
# rm ${trimmed}
# rm ${intersection}
# rm ${abc}
