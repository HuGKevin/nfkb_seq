#!/bin/bash

# MINT_MASTER.sh
#
# The main MintChIP  processing scripts for NFKB MintChIP project. 
#
# As currently configured, this script runs on Farnam.
#
# The following steps are performed:
#   1. Check that all required components are present
#   2. QC of raw reads with FastQC
#   3. Adapter trimming with pyadapter_trim.py
#   4. QC of trimmed reads with FastQC
#   5. Alignment of trimmed reads with Bowtie2
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

############ TRIM ADAPTERS #############
#!/bin/bash

#SBATCH --partition=general
#SBATCH --job-name=trim_mint%a
#SBATCH --cpus-per-task=1 --mem=15gb
#SBATCH -o /home/kh593/scratch60/nfkb_seq/logs/trim_mint%a.out
#SBATCH -e /home/kh593/scratch60/nfkb_seq/logs/trim_mint%a.err
#SBATCH --array=2

# Clear out environment of node and load conda environment
module purge
module load miniconda
source activate atac_env
echo "conda environment activated"
python --version

# File with metadata job array
index_file=/home/kh593/project/nfkb_seq/data/mint_copa.tsv

# Extract relevant arguments from the table
donor=$(awk -F'\t' -v row=${SLURM_ARRAY_TASK_ID} -v num=1 'FNR == row {print $num}' $index_file)
expt=$(awk -F'\t' -v row=${SLURM_ARRAY_TASK_ID} -v num=2 'FNR == row {print $num}' $index_file)
stim=$(awk -F'\t' -v row=${SLURM_ARRAY_TASK_ID} -v num=3 'FNR == row {print $num}' $index_file)
dnalib=$(awk -F'\t' -v row=${SLURM_ARRAY_TASK_ID} -v num=4 'FNR == row {print $num}' $index_file)
copa=$(awk -F'\t' -v row=${SLURM_ARRAY_TASK_ID} -v num=5 'FNR == row {print $num}' $index_file)

# Relevant directories
data_dir="/home/kh593/scratch60/nfkb_data/MintChIP"
src_dir="/home/kh593/project/nfkb_seq/src"
trimmed_reads_dir="/home/kh593/scratch60/nfkb_seq/trimmed_reads/mint"

# Intermediate files
fastq1_trim="${trimmed_reads_dir}/${copa}_R1.trim.fastq.gz"
fastq2_trim="${trimmed_reads_dir}/${copa}_R2.trim.fastq.gz"
 
# Trim adapter sequences:
python ${src_dir}/pyadapter_trim.py -a ${data_dir}/${copa}_R1.fastq.gz -b ${data_dir}/${copa}_R2.fastq.gz

mv ${copa}_R1.trim.fastq.gz ${trimmed_reads_dir}
mv ${copa}_R2.trim.fastq.gz ${trimmed_reads_dir}

echo "Adapters trimmed for ${copa}"

echo "Modules listed below for reference"

module list
#########################

### Below FastQC is optional
# Perform QC testing on raw reads with FastQC:
fastqc ${data_dir}/${R1} ${data_dir}/${R2} --outdir=$pretrim_qc_dir/atac --threads=16

echo 'Raw QC Completed'

# Perform QC testing on trimmed reads with FastQC:       
fastqc $fastq1_trim $fastq2_trim --outdir=$posttrim_qc_dir
echo 'Post-trim QC completed'


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
Rscript /home/kh593/project/nfkb_seq/src/peak_matrix.R


