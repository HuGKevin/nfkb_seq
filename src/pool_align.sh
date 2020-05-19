#!/bin/bash

#SBATCH --partition=scavenge
#SBATCH --job-name=pool_align%a
#SBATCH --cpus-per-task=16 --mem=45gb
#SBATCH -o /home/kh593/scratch60/nfkb_seq/logs/pool_align%a.out
#SBATCH -e /home/kh593/scratch60/nfkb_seq/logs/pool_align%a.err
#SBATCH --array=106,112,113,114,117,119,121,122,124,125,129,131,132,135,137,141,142,143,145,150,154,157,167,27,28,90

# Clear out environment of node and load conda environment
module purge

# Load additional necessary packages
module load SAMtools
# module load Biopython
# module load Python
module load Bowtie2
module load picard

python --version

# File with metadata job array
index_file="/home/kh593/project/nfkb_seq/data/atac_libs.tsv"
bt2_index="/gpfs/ysm/project/cotsapas/kh593/genomes/hg38/Bowtie2_index/hg38_index"
genome_seq="/gpfs/ysm/project/cotsapas/kh593/genomes/hg38/hg38.fa"

# Extract relevant arguments from the table
donor=$(awk -F'\t' -v row=${SLURM_ARRAY_TASK_ID} -v num=1 'FNR == row {print $num}' $index_file)
expt=$(awk -F'\t' -v row=${SLURM_ARRAY_TASK_ID} -v num=2 'FNR == row {print $num}' $index_file)
stim=$(awk -F'\t' -v row=${SLURM_ARRAY_TASK_ID} -v num=3 'FNR == row {print $num}' $index_file)
lib=$(awk -F'\t' -v row=${SLURM_ARRAY_TASK_ID} -v num=4 'FNR == row {print $num}' $index_file)

# Relevant directories
scratch_dir="/home/kh593/scratch60/nfkb_seq"
alignment_dir="${scratch_dir}/aligned_reads"
logs_dir="${scratch_dir}/logs"
alignment_stats_dir="${scratch_dir}/alignment_stats"

# Intermediate files
r1_pooled="/home/kh593/scratch60/nfkb_seq/pooled_reads/${lib}_R1_pooled.fastq.gz"
r2_pooled="/home/kh593/scratch60/nfkb_seq/pooled_reads/${lib}_R2_pooled.fastq.gz"
aligned_file="${alignment_dir}/${lib}.bam"
alignment_stats_file="${alignment_stats_dir}/${lib}.alignment.stats.txt"
duplicate_log="${logs_dir}/${lib}.atac.duplicates.log"
filtered_file="${alignment_dir}/${lib}.filtered.bam"
nodup_file="${alignment_dir}/${lib}.nodup.bam"
shifted_file="${alignment_dir}/${lib}.shifted.bam"
final_file="${alignment_dir}/${lib}.final.bam"
final_alignment_stats="${alignment_stats_dir}/${lib}.final.alignment.stats.txt"

# Pool R1s and R2s. 
r1s=$(grep "${lib}" /home/kh593/project/nfkb_seq/data/atac_copa.tsv |
	  cut -f5 |
	  sed "s/$/_R1.trim.fastq.gz/g" |
	  sed "s/^/\/home\/kh593\/scratch60\/nfkb_seq\/trimmed_reads\/atac\//g")
r2s=$(grep "${lib}" /home/kh593/project/nfkb_seq/data/atac_copa.tsv |
	  cut -f5 |
	  sed "s/$/_R2.trim.fastq.gz/g" |
	  sed "s/^/\/home\/kh593\/scratch60\/nfkb_seq\/trimmed_reads\/atac\//g")

cat ${r1s} > ${r1_pooled}
cat ${r2s} > ${r2_pooled}

echo "${lib} Pooled"

# Align reads with Bowtie2:
bowtie2 -p 16 -x $bt2_index -1 ${r1_pooled} -2 ${r2_pooled} |
    samtools sort -o $aligned_file -@ 16 -
samtools index $aligned_file
echo 'Reads aligned'

# Compile alignment statistics:
java -jar $EBROOTPICARD/picard.jar CollectAlignmentSummaryMetrics R=$genome_seq I=$aligned_file O=$alignment_stats_file
echo 'Alignment stats compiled'

# Filter reads for quality, both ends properly mapped, and from principal assembly:
samtools view -q 30 -f 0x2 -h -@ 16 $aligned_file chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY |
    samtools sort -o $filtered_file -@ 16 -
samtools index $filtered_file

# Remove precursor files
if [ -f $filtered_file ]; then
    rm ${aligned_file}
    rm ${aligned_file}.bai
fi

echo 'Reads filtered on quality, mapping, and source'

# Remove duplicate reads:
java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=$filtered_file O=$nodup_file M=$duplicate_log REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=SILENT
samtools index $nodup_file

# Remove precursor files
if [ -f $nodup_file ]; then
    rm ${filtered_file}
    rm ${filtered_file}.bai
fi

echo 'Duplicate reads removed'

# Call shift_reads.py to shift read cut sites:
python /home/kh593/project/nfkb_seq/src/shift_reads.py -o $shifted_file $nodup_file
samtools sort -o $final_file -@ 16 $shifted_file
samtools index $final_file

# Remove precursor files
if [ -f ${final_file} ]; then
    rm ${nodup_file}
    rm ${nodup_file}.bai
    rm ${shifted_file}
fi
echo 'Read cut sites shifted'

# Compile alignment statistics on shifted reads:
java -jar $EBROOTPICARD/picard.jar CollectAlignmentSummaryMetrics R=$genome_seq I=$final_file O=$final_alignment_stats
echo 'Shifted reads alignment stats compiled'
