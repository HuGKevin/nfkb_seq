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
#SBATCH --array=2-175

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

################# POOL AND ALIGN ##############
#!/bin/bash

#SBATCH --partition=general
#SBATCH --job-name=pool_align_mint%a
#SBATCH --cpus-per-task=16 --mem=20gb
#SBATCH -o /home/kh593/scratch60/nfkb_seq/logs/pool_align_mint%a.out
#SBATCH -e /home/kh593/scratch60/nfkb_seq/logs/pool_align_mint%a.err
#SBATCH --array=2-175

echo "Job ID: ${SLURM_JOB_ID}"
echo "Array ID: ${SLURM_JOB_ARRAY_ID}"

# Clear out environment of node and load conda environment
module purge

# Load additional necessary packages
module load SAMtools
module load Bowtie2
module load picard

python --version

# File with metadata job array
index_file="/home/kh593/project/nfkb_seq/data/mint_libs.tsv"
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
filtered_file="${alignment_dir}/${lib}.filtered.bam"

# if [ -f ${final_file} ]
# then
#     echo "Alignment already completed" 
#     exit 4
# fi

# Pool R1s and R2s. 
r1s=$(grep "${lib}" /home/kh593/project/nfkb_seq/data/mint_copa.tsv |
	  cut -f5 |
	  sed "s/$/_R1.trim.fastq.gz/g" |
	  sed "s/^/\/home\/kh593\/scratch60\/nfkb_seq\/trimmed_reads\/mint\//g")
r2s=$(grep "${lib}" /home/kh593/project/nfkb_seq/data/mint_copa.tsv |
	  cut -f5 |
	  sed "s/$/_R2.trim.fastq.gz/g" |
	  sed "s/^/\/home\/kh593\/scratch60\/nfkb_seq\/trimmed_reads\/mint\//g")

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
samtools view -q 30 -f 0x2 -h -@ 16 $aligned_file \
	 chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
	 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 \
	 chr21 chr22 chrX chrY |
    samtools sort -o $filtered_file -@ 16 -
samtools index $filtered_file

# Remove precursor files
if [ -f $filtered_file ]; then
    rm ${aligned_file}
    rm ${aligned_file}.bai
fi

echo 'Reads filtered on quality, mapping, and source'
######################################################

################# DEDUP AND STATS ##############
#!/bin/bash

#SBATCH --partition=general
#SBATCH --job-name=dedup_stats_mint%a
#SBATCH --cpus-per-task=1 --mem=50gb
#SBATCH -o /home/kh593/scratch60/nfkb_seq/logs/dedup_stats_mint%a.out
#SBATCH -e /home/kh593/scratch60/nfkb_seq/logs/dedup_stats_mint%a.err
#SBATCH --array=2-174

echo "Job ID: ${SLURM_JOB_ID}"
echo "Array ID: ${SLURM_JOB_ARRAY_ID}"

# Clear out environment of node and load conda environment
module purge

# Load additional necessary packages
module load SAMtools
module load Bowtie2
module load picard

python --version

# File with metadata job array
index_file="/home/kh593/project/nfkb_seq/data/mint_libs.tsv"
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
duplicate_log="${logs_dir}/${lib}.mint.duplicates.log"
filtered_file="${alignment_dir}/${lib}.filtered.bam"
final_file="${alignment_dir}/${lib}.final.bam"
final_alignment_stats="${alignment_stats_dir}/${lib}.final.alignment.stats.txt"

# Remove duplicate reads:
java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=$filtered_file O=$final_file M=$duplicate_log REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=SILENT
samtools index $final_file

# # Remove precursor files
# if [ -f $final_file ]; then
#     rm ${filtered_file}
#     rm ${filtered_file}.bai
# fi

echo 'Duplicate reads removed'

# Compile alignment statistics on shifted reads:
java -jar $EBROOTPICARD/picard.jar CollectAlignmentSummaryMetrics R=$genome_seq I=$final_file O=$final_alignment_stats
echo 'Alignment stats compiled'
##################################################################


################################################################################
###########################   Section 3: Find Peaks   ##########################
################################################################################

################## Downsample and call peaks #################
#!/bin/bash

#SBATCH --partition=general
#SBATCH --job-name=peakcall_mint%a
#SBATCH --cpus-per-task=8 --mem=15gb
#SBATCH -o /home/kh593/scratch60/nfkb_seq/logs/peakcall_mint%a.out
#SBATCH -e /home/kh593/scratch60/nfkb_seq/logs/peakcall_mint%a.err
#SBATCH --array=2-174

echo "Job ID: ${SLURM_JOB_ID}"
echo "Array ID: ${SLURM_JOB_ARRAY_ID}"

module purge
module load SAMtools
module load MACS2

index_file="/home/kh593/project/nfkb_seq/data/mint_libs.tsv"

### job information
donor=$(awk -F'\t' -v row=${SLURM_ARRAY_TASK_ID} -v num=1 'FNR == row {print $num}' $index_file)
expt=$(awk -F'\t' -v row=${SLURM_ARRAY_TASK_ID} -v num=2 'FNR == row {print $num}' $index_file)
stim=$(awk -F'\t' -v row=${SLURM_ARRAY_TASK_ID} -v num=3 'FNR == row {print $num}' $index_file)
lib=$(awk -F'\t' -v row=${SLURM_ARRAY_TASK_ID} -v num=4 'FNR == row {print $num}' $index_file)

# Intermediate files and directories
scratch_dir="/home/kh593/scratch60/nfkb_seq"
reads_dir="${scratch_dir}/aligned_reads"
peaks_dir="${scratch_dir}/peaks"
fullset="${scratch_dir}/${lib}.final.bam"
downsample="${scratch_dir}/${lib}.down.bam"
npeakfile="${peaks_dir}/${lib}.narrowPeaks.gz"
bpeakfile="${peaks_dir}/${lib}.broadPeaks.gz"
gpeakfile="${peaks_dir}/${lib}.gappedPeaks.gz"
FE_bedgraph="${lib}_FE.bdg"
Pval_bedgraph="${lib}_ppois.bdg"
FE_bigwig="${peaks_dir}/${lib}_FE.bigwig"
Pval_bigwig="${peaks_dir}/${lib}_ppois.bigwig"

# Necessary variables and files
target_depth=30000000
pval_thresh="0.01"
genome_sizes="/home/kh593/project/genomes/hg38/hg38_principal.chrom.sizes"

# compute total reads
total_reads=$(samtools view -@ 8 -c ${fullset})
echo "Target depth: ${target_depth}"
echo "Fullset size: ${total_reads}"

# compute fraction of reads given an input read depth
frac=$(awk -v down=$target_depth -v full=$total_reads 'BEGIN {frac=down/full;
if (frac > 1) {print 1} else {print frac}}')	  

# samtools view to downsample
if [ $frac -eq 1 ]
then
    echo "${lib} doesn't exceed downsample threshold"
    exit 25
else
    echo "Downsample fraction: ${frac}"
    samtools view -@ 8 -bs $frac $fullset > $downsample
    samtools index $downsample
fi

# Peak call atac-seq peaks using narrowpeaks from MACS2
### Call narrow peaks
macs2 callpeak \
      -t $downsample -f BAMPE -n $lib -p $pval_thresh \
      --nomodel -B --SPMR --keep-dup all --call-summits

### Call broad peaks
macs2 callpeak \
      -t $downsample -f BAMPE -n $lib -p $pval_thresh \
      --nomodel --broad --keep-dup all

### Sort peaks and compress them
## narrow
sort -k 8gr,8gr ${lib}_peaks.narrowPeak | \
    awk -v id=$lib 'BEGIN{OFS = "\t"}{$4 = id"_peak"NR ; print $0}' | \
    gzip -c > $npeakfile

## broad
sort -k 8gr,8gr ${lib}_peaks.broadPeak | \
    awk -v id=$lib 'BEGIN{OFS = "\t"}{$4 = id"_peak"NR ; print $0}' | \
    gzip -c > $bpeakfile

## gapped
sort -k 8gr,8gr ${lib}_peaks.gappedPeak | \
    awk -v id=$lib 'BEGIN{OFS = "\t"}{$4 = id"_peak"NR ; print $0}' | \
    gzip -c > $gpeakfile

### Convert peaks to signal tracks
## Fold enrichment track
macs2 bdgcmp -t ${lib}_treat_pileup.bdg \
      -c ${lib}_control_lambda.bdg \
      --o-prefix ${lib} -m FE

## Poisson pval track
macs2 bdgcmp -t ${lib}_treat_pileup.bdg \
      -c ${lib}_control_lambda.bdg \
      --o-prefix ${lib} -m ppois

## Sort bedgraphs
bedSort ${FE_bedgraph} ${FE_bedgraph}
bedSort ${Pval_bedgraph} ${Pval_bedgraph}

### Convert signal tracks to bigwigs
bedGraphToBigWig ${FE_bedgraph} $genome_sizes $FE_bigwig
bedGraphToBigWig ${Pval_bedgraph} $genome_sizes $Pval_bigwig

if [ -f "${FE_bigwig}" ] 
then
    rm -f ${FE_bedgraph}
fi

if [ -f "${Pval_bigwig}" ]
then
    rm -f ${Pval_bedgraph}
fi

### Remove unnecessary files
rm -f ${lib}_peaks.narrowPeak
rm -f ${lib}_peaks.broadPeak
rm -f ${lib}_peaks.gappedPeak
rm -f ${lib}_peaks.xls
rm -f ${lib}_treat_pileup.bdg
rm -f ${lib}_control_lambda.bdg

### Move summits to peak dir
mv ${lib}_summits.bed ${peaks_dir}
##############################################################

### Run MCL on each set of experiments
mcl_dir="${scratch_dir}/mcl"
mint_dir="${mcl_dir}/mint"

gunzip /home/kh593/scratch60/nfkb_seq/peaks/mint/*.gz

cat /home/kh593/scratch60/nfkb_seq/peaks/mint/*.broadPeaks > ${mint_dir}/mcl_compiled.bed

bedSort ${mint_dir}/mcl_compiled.bed ${mint_dir}/mcl_compiled.bed
awk '{print $0 >> $1".bed"}' ${mint_dir}/mcl_compiled.bed
mv chr*.bed ${mint_dir}

sbatch /home/kh593/project/nfkb_seq/src/mcl.sh

### Bind cluster output from MCL
sbatch /home/kh593/project/nfkb_seq/src/bind_clusters.sh

### Generate binding profile for each DNA library
Rscript /home/kh593/project/nfkb_seq/src/peak_matrix.R


