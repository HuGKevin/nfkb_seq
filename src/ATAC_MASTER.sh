#!/bin/bash

# ATAC_MASTER.sh
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
#   5. Alignment of trimmed reads with Bowtie2
#   6. Compilation of alignment summary statistics by Picard tools
#   7. Filtering to include only high quality singly mapped reads mapped to primary assembly
#   8. Filtering to remove optical duplicates
#   9. Shift read cut sites with shift_reads.py
#   10. Compilation of alignment summary statistics by Picard tools
#   11. Downsample each library to 30 million reads
#   12. Call peaks using the narrowPeak settings in MACS2
#   13. Process called peaks for parallel Markov Clustering
#   14. Run Markov Clustering
#   15. Convert MCL output to bed file and generate a presence/absence matrix

################################################################################
#############################   Section 1: Setup   #############################
################################################################################

# This section defines executable paths and locations of raw data. Any missing
# output directories are created.
# Most of this section exists as a legacy from Noah/Matt's code. It's not 100% necessary. 

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
peak_dir="${scratch_dir}/peaks"
mcl_dir="${scratch_dir}/mcl"

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

# Genome index:
bt2_index="/gpfs/ysm/project/cotsapas/kh593/genomes/hg38/Bowtie2_index/hg38_index"
genome_seq="/gpfs/ysm/project/cotsapas/kh593/genomes/hg38/hg38.fa"

# Check that all essential components are present:
essential_files=( ${bt2_index}.1.bt2 $genome_seq \
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

# This section contains job array scripts for QC, pre-alignment processing, and alignment

######### QC AND TRIM ADAPTERS #########
#!/bin/bash

#SBATCH --partition=general
#SBATCH --job-name=trim_atac%a
#SBATCH --cpus-per-task=1 --mem=35gb
#SBATCH -o /home/kh593/scratch60/nfkb_seq/logs/trim_atac%a.out
#SBATCH -e /home/kh593/scratch60/nfkb_seq/logs/trim_atac%a.err
#SBATCH --array=2-1145

# Set up environment
module purge
module load miniconda
source activate atac_env
echo "conda environment activated"
python --version

# File with metadata job array
index_file=/home/kh593/project/nfkb_seq/data/atac_copa.tsv

# Extract relevant arguments from the metadata
donor=$(awk -F'\t' -v row=${SLURM_ARRAY_TASK_ID} -v num=1 'FNR == row {print $num}' $index_file)
expt=$(awk -F'\t' -v row=${SLURM_ARRAY_TASK_ID} -v num=2 'FNR == row {print $num}' $index_file)
stim=$(awk -F'\t' -v row=${SLURM_ARRAY_TASK_ID} -v num=3 'FNR == row {print $num}' $index_file)
dnalib=$(awk -F'\t' -v row=${SLURM_ARRAY_TASK_ID} -v num=4 'FNR == row {print $num}' $index_file)
copa=$(awk -F'\t' -v row=${SLURM_ARRAY_TASK_ID} -v num=5 'FNR == row {print $num}' $index_file)

# Relevant directories
data_dir="/home/kh593/scratch60/nfkb_data/ATACseq"
src_dir="/home/kh593/project/nfkb_seq/src"
trimmed_reads_dir="/home/kh593/scratch60/nfkb_seq/trimmed_reads/atac"

# Trim adapter sequences:
python ${src_dir}/pyadapter_trim.py -a ${data_dir}/${copa}_R1.fastq.gz -b ${data_dir}/${copa}_R2.fastq.gz

mv ${copa}_R1.trim.fastq.gz ${trimmed_reads_dir}
mv ${copa}_R2.trim.fastq.gz ${trimmed_reads_dir}

echo "Adapters trimmed for ${copa}"

echo "Modules listed below for reference"

module list
#########################

### Below FastQC is optional
# Perform QC testing on raw reads with FastQC (multithreaded to 16 cores):
fastqc ${data_dir}/${R1} ${data_dir}/${R2} --outdir=$pretrim_qc_dir/atac --threads=16

echo 'Raw QC Completed'

# Perform QC testing on trimmed reads with FastQC:       
fastqc $fastq1_trim $fastq2_trim --outdir=$posttrim_qc_dir
echo 'Post-trim QC completed'


######## POOL AND ALIGN ############
#!/bin/bash

#SBATCH --partition=general
#SBATCH --job-name=pool_align_atac%a
#SBATCH --cpus-per-task=16 --mem=20gb
#SBATCH -o /home/kh593/scratch60/nfkb_seq/logs/pool_align_atac%a.out
#SBATCH -e /home/kh593/scratch60/nfkb_seq/logs/pool_align_atac%a.err
#SBATCH --array=2-183

# Set up environment
module purge
module load SAMtools
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
filtered_file="${alignment_dir}/${lib}.filtered.bam"

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

# Align reads with Bowtie2
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
#############################################

######## DEDUP AND STATS  ############
#!/bin/bash

#SBATCH --partition=general
#SBATCH --job-name=dedup_stats_atac%a
#SBATCH --cpus-per-task=8 --mem=30gb
#SBATCH -o /home/kh593/scratch60/nfkb_seq/logs/dedup_stats_atac%a.out
#SBATCH -e /home/kh593/scratch60/nfkb_seq/logs/dedup_stats_atac%a.err
#SBATCH --array=2-183

# Set up environment
module purge
module load SAMtools
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
duplicate_log="${logs_dir}/${lib}.atac.duplicates.log"
filtered_file="${alignment_dir}/${lib}.filtered.bam"
nodup_file="${alignment_dir}/${lib}.nodup.bam"
shifted_file="${alignment_dir}/${lib}.shifted.bam"
final_file="${alignment_dir}/${lib}.final.bam"
final_alignment_stats="${alignment_stats_dir}/${lib}.final.alignment.stats.txt"

echo $lib

# Remove duplicate reads:
java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=$filtered_file O=$nodup_file M=$duplicate_log REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=SILENT
samtools index $nodup_file

echo 'Duplicate reads removed'

# Call shift_reads.py to shift read cut sites:
python /home/kh593/project/nfkb_seq/src/shift_reads2.py -o $shifted_file $nodup_file
samtools sort -o $final_file -@ 4 $shifted_file
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
###########################################

################################################################################
###########################   Section 3: Call Peaks   ##########################
################################################################################

################## Downsample and call peaks #################
#!/bin/bash

#SBATCH --partition=general
#SBATCH --job-name=peakcall_atac%a
#SBATCH --cpus-per-task=8 --mem=15gb
#SBATCH -o /home/kh593/scratch60/nfkb_seq/logs/peakcall_atac%a.out
#SBATCH -e /home/kh593/scratch60/nfkb_seq/logs/peakcall_atac%a.err
#SBATCH --array=2-183

echo "Job ID: ${SLURM_JOB_ID}"
echo "Array ID: ${SLURM_ARRAY_TASK_ID}"

# Set up environment
module purge
module load SAMtools
module load MACS2

index_file="/home/kh593/project/nfkb_seq/data/atac_libs.tsv"

### job information
donor=$(awk -F'\t' -v row=${SLURM_ARRAY_TASK_ID} -v num=1 'FNR == row {print $num}' $index_file)
expt=$(awk -F'\t' -v row=${SLURM_ARRAY_TASK_ID} -v num=2 'FNR == row {print $num}' $index_file)
stim=$(awk -F'\t' -v row=${SLURM_ARRAY_TASK_ID} -v num=3 'FNR == row {print $num}' $index_file)
lib=$(awk -F'\t' -v row=${SLURM_ARRAY_TASK_ID} -v num=4 'FNR == row {print $num}' $index_file)

# Intermediate files and directories
scratch_dir="/home/kh593/scratch60/nfkb_seq"
reads_dir="${scratch_dir}/aligned_reads"
peaks_dir="${scratch_dir}/peaks"
fullset="${reads_dir}/${lib}.final.bam"
downsample="${reads_dir}/${lib}.down.bam"
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

# Compute total reads
total_reads=$(samtools view -@ 8 -c ${fullset})
echo "Target depth: ${target_depth}"
echo "Fullset size: ${total_reads}"

# Compute fraction of reads given an input read depth
frac=$(awk -v down=$target_depth -v full=$total_reads 'BEGIN {frac=down/full;
if (frac > 1) {print 1} else {print frac}}')	  

# Downsample to depth if library is sufficiently deep
if [ $frac -eq 1 ]
then
    echo "${lib} doesn't exceed downsample threshold"
    exit 25
else
    echo "Downsample fraction: ${frac}"
    samtools view -@ 8 -bs $frac $fullset > $downsample
    samtools index $downsample
fi

# Peak call atac-seq peaks MACS2
## Call narrow peaks
macs2 callpeak \
      -t $downsample -f BAMPE -n $lib -p $pval_thresh \
      --nomodel -B --SPMR --keep-dup all --call-summits

## Call broad peaks
macs2 callpeak \
      -t $downsample -f BAMPE -n $lib -p $pval_thresh \
      --nomodel --broad --keep-dup all

## Sort peaks and compress them
### Narrow
sort -k 8gr,8gr ${lib}_peaks.narrowPeak | \
    awk -v id=$lib 'BEGIN{OFS = "\t"}{$4 = id"_peak"NR ; print $0}' | \
    gzip -c > $npeakfile

### Broad
sort -k 8gr,8gr ${lib}_peaks.broadPeak | \
    awk -v id=$lib 'BEGIN{OFS = "\t"}{$4 = id"_peak"NR ; print $0}' | \
    gzip -c > $bpeakfile

### Gapped
sort -k 8gr,8gr ${lib}_peaks.gappedPeak | \
    awk -v id=$lib 'BEGIN{OFS = "\t"}{$4 = id"_peak"NR ; print $0}' | \
    gzip -c > $gpeakfile

# Convert peaks to signal tracks
## Fold enrichment track
macs2 bdgcmp -t ${lib}_treat_pileup.bdg \
      -c ${lib}_control_lambda.bdg \
      --o-prefix ${lib} -m FE

## Poisson pval track
macs2 bdgcmp -t ${lib}_treat_pileup.bdg \
      -c ${lib}_control_lambda.bdg \
      --o-prefix ${lib} -m ppois

# Sort bedgraphs
bedSort ${FE_bedgraph} ${FE_bedgraph}
bedSort ${Pval_bedgraph} ${Pval_bedgraph}

# Convert signal tracks to bigwigs
bedGraphToBigWig ${FE_bedgraph} $genome_sizes $FE_bigwig
bedGraphToBigWig ${Pval_bedgraph} $genome_sizes $Pval_bigwig

# Remove precursor files if bigwig conversion successful
if [ -f "${FE_bigwig}" ] 
then
    rm -f ${FE_bedgraph}
fi

if [ -f "${Pval_bigwig}" ]
then
    rm -f ${Pval_bedgraph}
fi

# Remove unnecessary files
rm -f ${lib}_peaks.narrowPeak
rm -f ${lib}_peaks.broadPeak
rm -f ${lib}_peaks.gappedPeak
rm -f ${lib}_peaks.xls
rm -f ${lib}_treat_pileup.bdg
rm -f ${lib}_control_lambda.bdg

# Move summits to peak dir
mv ${lib}_summits.bed ${peaks_dir}

#####################################################

################################################################################
#########################   Section 4: Cluster Peaks  ##########################
################################################################################

##################### Prepare files for MCL #######################
mcl_dir="${scratch_dir}/mcl"
atac_dir="${mcl_dir}/atac"

# Combine all peaks into a single file
gunzip /home/kh593/scratch60/nfkb_seq/peaks/atac/*.gz 
echo -n > ${atac_dir}/mcl_compiled.bed 

## Exclude tech-dev donors from clustering
while read donor expt stim lib
do
    if [ $donor == "donor" ] ||
	   [ $donor == "TB5728" ] || 
	   [ $donor == "TB0611" ] ||
	   [ $donor == "TB6578" ]
    then
	echo "skipped ${lib}"
	continue
    else
	cat ${peak_dir}/atac/${lib}.narrowPeaks >> ${atac_dir}/mcl_compiled.bed
    fi
    
done < ~/project/nfkb_seq/data/atac_libs.tsv

# Split filtered peaks by chromosome to make MCL parallel.
bedSort ${atac_dir}/mcl_compiled.bed ${atac_dir}/mcl_compiled.bed 
awk '{print $0 >> $1".bed"}' ${atac_dir}/mcl_compiled.bed 
mv chr*.bed ${atac_dir}

########################################################################

################### Run MCL on each chromosome #####################
#!/bin/bash

#SBATCH --partition=general
#SBATCH --job-name=mcl_atac%a
#SBATCH --cpus-per-task=10 --mem=64gb
#SBATCH -o /home/kh593/scratch60/nfkb_seq/logs/mcl_atac%a.out
#SBATCH -e /home/kh593/scratch60/nfkb_seq/logs/mcl_atac%a.err
#SBATCH --array=1-24

# Set up environment
module purge
module load MCL
module load BEDTools

if [ $SLURM_ARRAY_TASK_ID -eq 23 ]
then
    chr=X
elif [ $SLURM_ARRAY_TASK_ID -eq 24 ]
then
    chr=Y
else
    chr=$SLURM_ARRAY_TASK_ID
fi

mcl_dir="/home/kh593/scratch60/nfkb_seq/mcl"
file="${mcl_dir}/atac/chr${chr}.bed"
trimmed="${mcl_dir}/atac/chr${chr}_trimmed.bed"
intersection="${mcl_dir}/atac/chr${chr}_int.bed"
abc="${mcl_dir}/atac/chr${chr}_input.abc"
final="${mcl_dir}/atac/chr${chr}_clusters.txt"

# Cut out excess columns
cut -f1,2,3,4 $file > ${trimmed}

# Intersect against self and report overlaps and size of overlap
intersectBed -a ${trimmed} -b ${trimmed} -wo -sorted > ${intersection}

# Remove unnecessary columns
cut -f4,8,9 ${intersection} > ${abc}

# Run MCL (for future note, can use -te for multithreading)
mcl ${abc} --abc -o ${final} -te 10 

# Remove intermediate files
rm ${trimmed}
rm ${intersection}
rm ${abc}

########################################################################

############### Convert MCL output to BED file of clusters #################
#!/bin/bash
# Usage: sbatch bind_clusters.sh

#SBATCH --job-name=bind_atac_clusters%a
#SBATCH --partition=general
#SBATCH --mem=25gb --cpus-per-task=1
#SBATCH -e /home/kh593/scratch60/nfkb_seq/logs/bind_clusters%a.err
#SBATCH -o /home/kh593/scratch60/nfkb_seq/logs/bind_clusters%a.out
#SBATCH --array=1-24

# Set up environment
module purge
module load R/3.6.1-foss-2018b

index="/home/kh593/project/nfkb_seq/data/mcl_index.tsv"

chr=$(awk -F"\t" -v row=${SLURM_ARRAY_TASK_ID} -v num=1 'FNR == row {print $num}' $index)
expt=$(awk -F"\t" -v row=${SLURM_ARRAY_TASK_ID} -v num=2 'FNR == row {print $num}' $index)

Rscript /home/kh593/project/nfkb_seq/src/bind_clusters.R atac ${chr}

echo "Job completed"

########################################################################

################################################################################
########################   Section 5: Cluster Analysis  ########################
################################################################################

# Relevant libraries
module load BEDTools/2.29.2-foss-2018b
module load MCL/14.137-foss-2016b

# Relevant directories
project_dir="/home/kh593/project/nfkb_seq"
scratch_dir="/home/kh593/scratch60/nfkb_seq"
peak_dir="${scratch_dir}/peaks"
mcl_dir="${scratch_dir}/mcl"
analysis_dir="${scratch_dir}/analysis"
overlap_dir="${analysis_dir}/overlaps"

# Relevant input files
atac_clusters="${mcl_dir}/final_atac.bed"
mint_clusters="${mcl_dir}/final_mint.bed"
index="${project_dir}/data/called_libs.tsv"

# Relevant output files
consensus_rep_output="${analysis_dir}/consensus_replications.txt"

# Compare individual library peaksets
echo -e "donor\texpt\tstim\tlib\tpeaks_reps" > ${consensus_rep_output}

while read donor expt stim lib
do
    if [ $donor == "Donor" ] || [ $expt == "WCE" ]
    then
	continue
    fi

    if [ $expt == "H3K27ac" ]
    then
	cluster=${mint_clusters}
	peaks="${peak_dir}/mint/${lib}.broadPeaks"
    elif [ $expt == "ATAC" ]
    then
	cluster=${atac_clusters}
	peaks="${peak_dir}/atac/${lib}.narrowPeaks"
    fi
        
    ### Number of peaks replicated
    intersections=$(intersectBed -a $cluster -b $peaks -u | wc | awk {'print $1'})
    echo -e "${donor}\t${expt}\t${stim}\t${lib}\t${intersections}" >> ${consensus_rep_output} ## This is actually just a measure of peak density in the clusters, since every peak is represented in a cluster...

    ### Also create a bed file that shows how many peaks are in each cluster in each library
    intersectBed -a $cluster -b $peaks -c > ${overlap_dir}/${lib}_count.bed # Last column is number of peaks from this lib in that cluster

    ### Overlaps of peak replications
    if [ $expt == "ATAC" ]
    then
	intersectBed -a $cluster -b $peaks -wo | cut -f1,2,3,7,8,9,17 | sort | uniq > ${overlap_dir}/${lib}_overlaps.bed
    elif [ $expt == "H3K27ac" ]
    then
	intersectBed -a $cluster -b $peaks -wo | cut -f1,2,3,7,8,9,16 | sort | uniq > ${overlap_dir}/${lib}_overlaps.bed
    fi
	
    bedSort ${overlap_dir}/${lib}_overlaps.bed ${overlap_dir}/${lib}_overlaps.bed

    echo "${donor} ${expt} ${stim} ${lib}"

done < ${index}

# Generate presence/absence matrix on MCL clusters
Rscript /home/kh593/project/nfkb_seq/src/peak_matrix.R

########### Filter out clusters below a certain class ##############


