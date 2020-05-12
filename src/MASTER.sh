#!/bin/bash

# MASTER.sh
#
# The main ATAC-seq processing scripts for NFKB ATAC-seq project. 
#
# As currently configured, this script runs on Farnam.
#
# The following steps are performed:
#   1. Check that all required components are present
#   2. QC of raw reads with FastQC
#   3. Adapter trimming with pyadapter_trim.py
#   4. QC of trimmed reads with FastQC
#   5. Alignment of trimmed reads with HISAT2
#   6. Compilation of alignment summary statistics by Picard tools
#   7. Filtering to include only singly mapped reads mapped to primary assembly
#   8. Filtering to remove optical duplicates
#   9. Shift read cut sites with shift_reads.py
#   10. Compilation of alignment summary statistics by Picard tools


################################################################################
#############################   Section 1: Setup   #############################
################################################################################

# This section defines executable paths and locations of raw data. Any missing
# output directories are created.

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
	echo "Creating directory ${directory}"
	mkdir -p $directory
    fi
done

# The sample table:
sample_table="${base_dir}/data/file_index.csv"

# Genome index:
bt2_index="/gpfs/ysm/project/cotsapas/kh593/genomes/hg38/Bowtie2_index/hg38_index"
genome_seq="/gpfs/ysm/project/cotsapas/kh593/genomes/hg38/hg38.fa"

# Check that all essential components are present:
essential_files=( $sample_table ${bt2_index}.1.bt2 $genome_seq \
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


################################################################################
##########################   Section 2: Align Reads   ##########################
################################################################################

# This section submits multiple slurm scripts for QC and alignment.

# Get length of file index (assuming it has header)
projsize=$(wc -l ${sample_table} | cut -f1 -d" ")

echo "Creating alignment job submission script."

  echo \
"#!/bin/bash

#SBATCH --partition=general
#SBATCH --job-name=align_reads%a
#SBATCH --cpus-per-task=16 --mem=32gb
#SBATCH -o ${logs_dir}/%a.out
#SBATCH -e ${logs_dir}/%a.err
#SBATCH --array=2-${projsize}

# Clear out environment of node and load conda environment
module purge
module load miniconda
source activate atac_env

# Load additional necessary packages
module load SAMtools
module load Biopython
module load Python
module load picard
module load FastQC
module load Bowtie2

# File with metadata job array
index_file=${sample_table}

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
echo 'Shifted reads alignment stats compiled'" > ${slurm_dir}/2_align_reads_array.sh

echo "Alignment job submission script created."
  
################################################################################
###########################   Section 3: Find Peaks   ##########################
################################################################################

# Directories we need
peaks_dir="${scratch_dir}/indiv_peaks"
consensus_dir="${scratch_dir}/cons_peaks"
pooled_dir="${scratch_dir}/pooled_reads"

# Create any missing directories:
directories=( $peaks_dir $consensus_dir )
for directory in ${directories[@]}; do
  if [ ! -d $directory ]; then
    echo "Creating directory ${directory}"
    mkdir -p $directory
  fi
done

#### Pool all the reads into their respective dna_libararies

sbatch /home/kh593/project/nfkb_seq/src/pool_reads.sh

### Downsample to 30000000 reads

sbatch /home/kh593/project/nfkb_seq/src/downsample.sh

### Do peak calling

sbatch /home/kh593/project/nfkb_seq/src/peak_call.sh

### Run MCL on each set of experiments
mcl_dir="${scratch_dir}/results/MCL"
atac_dir="${mcl_dir}/atac"
mint_dir="${mcl_dir}/mint"

gunzip /home/kh593/scratch60/nfkb_seq/results/peak_call/beds/atac/*.gz

for file in /home/kh593/scratch60/nfkb_seq/results/peak_call/beds/atac/*.narrowPeaks
do
    base=$(basename "${file}" .narrowPeaks)

    sed -i -e "s/Peak_\([0-9]*\)/Peak_${base}_\1/g" ${file}
done

cat /home/kh593/scratch60/nfkb_seq/results/peak_call/beds/atac/*.narrowPeaks > ${atac_dir}/mcl_compiled.bed

bedSort ${atac_dir}/mcl_compiled.bed ${atac_dir}/mcl_compiled.bed
awk '{print $0 >> $1".bed"}' ${atac_dir}/mcl_compiled.bed
mv chr*.bed ${atac_dir}

sbatch /home/kh593/project/nfkb_seq/src/mcl.sh

gunzip /home/kh593/scratch60/nfkb_seq/results/peak_call/beds/mint/*broadPeaks.gz

for file in /home/kh593/scratch60/nfkb_seq/results/peak_call/beds/mint/*.broadPeaks
do
    base=$(basename "${file}" .broadPeaks)

    sed -i -e "s/Peak_\([0-9]*\)/Peak_${base}_\1/g" ${file}
done

cat /home/kh593/scratch60/nfkb_seq/results/peak_call/beds/mint/*.broadPeaks > ${mint_dir}/mcl_compiled.bed

bedSort ${mint_dir}/mcl_compiled.bed ${mint_dir}/mcl_compiled.bed
awk '{print $0 >> $1".bed"}' ${mint_dir}/mcl_compiled.bed
mv chr*.bed ${mint_dir}

sbatch /home/kh593/project/nfkb_seq/src/mcl.sh

### Bind cluster output from MCL
sbatch /home/kh593/project/nfkb_seq/src/bind_clusters.sh

### Generate binding profile for each DNA library

