#!/bin/bash

#SBATCH --partition=general
#SBATCH --job-name=align_reads%a
#SBATCH --cpus-per-task=16 --mem=32gb
#SBATCH -o /home/kh593/scratch60/nfkb_seq/logs/%a.out
#SBATCH -e /home/kh593/scratch60/nfkb_seq/logs/%a.err
#SBATCH --array=2

################################################################################
#############################   Section 1: Setup   #############################
################################################################################

# This section defines executable paths and locations of raw data. Any missing
# output directories causes the job to stop and produce an error.

missing=0 # Stores whether any critical components are missing

# Raw data and working directories:
base_dir="/home/kh593/project/nfkb_seq"
src_dir="${base_dir}/src"
scratch_dir="/home/kh593/scratch60/nfkb_seq"
data_dir="/home/kh593/scratch60/nfkb_data"

# Output directories:
slurm_dir="${scratch_dir}/slurm"
logs_dir="${scratch_dir}/logs"
doc_dir="${scratch_dir}/doc"
results_dir="${base_dir}/results"

# Scratch directories:
pretrim_qc_dir="${scratch_dir}/pretrim_fastqc"
posttrim_qc_dir="${scratch_dir}/posttrim_fastqc"
trimmed_reads_dir="${scratch_dir}/trimmed_reads"
alignment_dir="${scratch_dir}/aligned_reads"
alignment_stats_dir="${scratch_dir}/alignment_stats"
qc_stats_dir="${scratch_dir}/qc_stats"

# Create any missing directories:
directories=( $base_dir $data_dir $src_dir $scratch_dir \
			$logs_dir $doc_dir $slurm_dir $pretrim_qc_dir \
			$posttrim_qc_dir $trimmed_reads_dir $alignment_dir $alignment_stats_dir \
			$qc_stats_dir $results_dir )
for directory in ${directories[@]}; do
    if [ ! -d $directory ]; then
	echo "Cannot find ${directory}"
	missing=1
    fi
done

# The sample table:
file_index="${base_dir}/data/file_index.csv"

# Genome index:
bt2_index="/gpfs/ysm/project/cotsapas/kh593/genomes/hg38/Bowtie2_index/hg38_index"
genome_seq="/gpfs/ysm/project/cotsapas/kh593/genomes/hg38/hg38.fa"

# Check that all essential components are present:
essential_files=( $file_index ${bt2_index}.1.bt2 $genome_seq \
                  ${src_dir}/pyadapter_trim.py ${src_dir}/shift_reads.py )
for file in ${essential_files[@]}; do
  if [ ! -f $file ]; then
    echo "Cannot find ${file}"
    missing=1
  fi
done

# Die if any required components missing:
if [ $missing = 1 ]; then
    echo "ERROR: Cannot find all required components."
    exit -1
fi

# Clear out environment of node and load conda environment
module purge
module load miniconda
source activate atac_env

# Load additional necessary packages
module load SAMtools
module load picard
module load FastQC
module load Bowtie2

# Extract relevant arguments from the table
fid=$(awk -F'\t' -v row=${SLURM_ARRAY_TASK_ID} -v num=1 'FNR == row {print $num}' $index_file)
R1=$(awk -F'\t' -v row=${SLURM_ARRAY_TASK_ID} -v num=2 'FNR == row {print $num}' $index_file)
R2=$(awk -F'\t' -v row=${SLURM_ARRAY_TASK_ID} -v num=3 'FNR == row {print $num}' $index_file)
expt=$(awk -F'\t' -v row=${SLURM_ARRAY_TASK_ID} -v num=4 'FNR == row {print $num}' $index_file)

# Intermediate files
fastq1_trim=${trimmed_reads_dir}/`echo $R1 | sed 's/.fastq.gz/.trim.fastq.gz/'`
fastq2_trim=${trimmed_reads_dir}/`echo $R2 | sed 's/.fastq.gz/.trim.fastq.gz/'`
aligned_file=${alignment_dir}/${fid}.bam
alignment_stats_file=${alignment_stats_dir}/${fid}.alignment.stats.txt
duplicate_log=${logs_dir}/${fid}.atac.duplicates.log
filtered_file=${alignment_dir}/${fid}.filtered.bam
nodup_file=${alignment_dir}/${fid}.nodup.bam
shifted_file=${alignment_dir}/${fid}.shifted.bam
final_file=${alignment_dir}/${fid}.final.bam
final_alignment_stats=${alignment_stats_dir}/${fid}.final.alignment.stats.txt
 
################################################################################
##########################   Section 2: Align Reads   ##########################
################################################################################

# This section runs QC on the raw reads, trims them, aligns them, and filters for quality.

# Perform QC testing on raw reads with FastQC:
fastqc ${data_dir}/${expt}/${R1} \\
       ${data_dir}/${expt}/${R2} \\
       --outdir=$pretrim_qc_dir \\
       --threads=16

echo 'Raw QC Completed'

# Trim adapter sequences:
python ${src_dir}/pyadapter_trim.py \\
-a ${data_dir}/${expt}/${R1} \\
-b ${data_dir}/${expt}/${R2}
       mv `echo $R1 | sed 's/.fastq.gz/.trim.fastq.gz/'` \\
       $trimmed_reads_dir
       mv `echo $R2 | sed 's/.fastq.gz/.trim.fastq.gz/'` \\
       $trimmed_reads_dir
       echo 'Adapters trimmed'

# Perform QC testing on trimmed reads with FastQC:       
fastqc $fastq1_trim \\
       $fastq2_trim \\
       --outdir=$posttrim_qc_dir
echo 'Post-trim QC completed'

# Align reads with Bowtie2:
bowtie2 -p 16 -x $bt2_index \\
       -1 $fastq1_trim \\
       -2 $fastq2_trim | \\
  samtools sort -o $aligned_file -@ 16 -
samtools index $aligned_file
echo 'Reads aligned'

# Compile alignment statistics:
java -jar \$EBROOTPICARD/picard.jar \\
          CollectAlignmentSummaryMetrics \\
          R=$genome_seq \\
          I=$aligned_file \\
          O=$alignment_stats_file
	  echo 'Alignment stats compiled'

# Filter reads for quality, both ends properly mapped, and from principal assembly:
samtools view -q 30 -f 0x2 -h -@ 16 \\
  $aligned_file \\
  chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 \\
  chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY | \\
  samtools sort -o $filtered_file -@ 16 -
samtools index $filtered_file

# Remove precursor files
if [ -f $sorted_aligned_file_filtered ]; then
  rm $sorted_aligned_file
  rm ${sorted_aligned_file}.bai
fi

echo 'Reads filtered on quality, mapping, and source'

# Remove duplicate reads:
java -jar \$EBROOTPICARD/picard.jar \\
          MarkDuplicates \\
          I=$filtered_file \\
          O=$nodup_file \\
          M=$duplicate_log \\
          REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=SILENT
samtools index $nodup_file

# Remove precursor files
if [ -f $sorted_aligned_file_nodup ]; then
  rm $sorted_aligned_file_filtered
  rm ${sorted_aligned_file_filtered}.bai
fi

echo 'Duplicate reads removed'

# Call shift_reads.py to shift read cut sites:
python ${src_dir}/shift_reads.py -o $shifted_file \\
  $nodup_file
samtools sort -o $final_file -@ 16 $shifted_file
samtools index $final_file

# Remove precursor files
if [ -f $sorted_aligned_file_shifted ]; then
  rm $sorted_aligned_file_nodup
  rm ${sorted_aligned_file_nodup}.bai
fi
echo 'Read cut sites shifted'

# Compile alignment statistics on shifted reads:
java -jar \$EBROOTPICARD/picard.jar \\
  CollectAlignmentSummaryMetrics \\
  R=$genome_seq \\
  I=$final_file \\
  O=$final_alignment_stats
echo 'Shifted reads alignment stats compiled'
