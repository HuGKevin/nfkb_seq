#!/bin/bash

#SBATCH --partition=general
#SBATCH --job-name=call_peaks.%a
#SBATCH --cpus-per-task=8
#SBATCH --mem=10gb
#SBATCH -o /home/kh593/scratch60/nfkb_seq/logs/call_peaks%a.out
#SBATCH -e /home/kh593/scratch60/nfkb_seq/logs/call_peaks%a.err
#SBATCH --array=1

module purge
module load miniconda
source activate atac_env

module load MACS2/2.2.7.1-foss-2018b-Python-3.7.0
module load BEDTools/2.29.2-foss-2018b

index_file="/home/kh593/project/nfkb_seq/data/peakcall_array.tsv"

### job information
donor=$(awk -F'\t' -v row=${SLURM_ARRAY_TASK_ID} -v num=1 'FNR == row {print $num}' $index_file)
expt=$(awk -F'\t' -v row=${SLURM_ARRAY_TASK_ID} -v num=2 'FNR == row {print $num}' $index_file)
stim=$(awk -F'\t' -v row=${SLURM_ARRAY_TASK_ID} -v num=3 'FNR == row {print $num}' $index_file)
lib=$(awk -F'\t' -v row=${SLURM_ARRAY_TASK_ID} -v num=4 'FNR == row {print $num}' $index_file)

### Relevant directories
scratch_dir="/home/kh593/nfkb_seq"
reads_dir="${scratch_dir}/downsampled"
results_dir="${scratch_dir}/results"
peakcall_dir="${results_dir}/peak_call"
bigwig_dir="${peakcall_dir}/bigwigs"
summits_dir="${peakcall_dir}/summits"
atac_dir="${peakcall_dir}/beds/atac"
mint_dir="${peakcall_dir}/beds/mint"

### Repeat parameters for logs
echo $donor
echo $expt
echo $stim
echo $lib

## Call peaks depending on experiment
if [ $expt == "ATAC" ]
then
    /home/kh593/project/nfkb_seq/src/call_peaks_atac.sh "${reads_dir}/${lib}.bam"
    
    mv ${lib}*narrowPeaks.gz ${atac_dir}
    mv ${lib}*.gz ${atac_dir}/otherpeaks
    mv ${lib}*.bigwig ${bigwig_dir}
    mv ${lib}_summits.bed ${summits_dir}
elif [ $expt == "H3K27ac" ] 
then
    /home/kh593/project/nfkb_seq/src/call_peaks_mint.sh "${reads_dir}/${lib}.bam"

    mv ${lib}*broadPeaks.gz ${mint_dir}
    mv ${lib}*gappedPeaks.gz ${mint_dir}/otherpeaks/
    mv ${lib}*.bigwig ${bigwig_dir}
elif [ $expt == "WCE" ]
then
    echo "${lib} is whole cell extract; dunno what to do about that right now"
else
    echo "Invalid experiment type"
    exit 999
fi

echo "Peak calling complete"

     
    
