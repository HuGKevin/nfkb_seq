#!/bin/bash
#SBATCH --partition=scavenge
#SBATCH --job-name=Compile.alignment.stats
#SBATCH --ntasks=1
#SBATCH --mem=5gb
#SBATCH -o /home/kh593/scratch60/nfkb_seq/logs/stat_compile.out
#SBATCH -e /home/kh593/scratch60/nfkb_seq/logs/stat_compile.err

# This script compiles qc metrics on all of our experiments into distinct files for each type of metric. 
# The overall structure of the script will be to unpack data, then to repack into qc-metric specific files, looping through all our samples for each metric. 

# Raw data and working directories:
base_dir="/home/kh593/project/nfkb_seq"
src_dir="${base_dir}/src"
scratch_dir="/home/kh593/scratch60/nfkb_seq"
data_dir="/home/kh593/scratch60/nfkb_data"
atac_dir="${data_dir}/atac"
mint_dir="${data_dir}/mint"

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

# final qc output directory:
finalqc_dir="${scratch_dir}/final_qc"

module load R

file_index="${base_dir}/data/file_index.csv"
file_index_header="Uid R1 R2 type copa"

while IFS=, read ${file_index_header}
do
  # Skip the header row and any samples that should be excluded:
  [ "${Uid}" == "UID" ] && continue

  echo "Processing sample ${copa}"
  
  unzip -p ${pretrim_qc_dir}/${copa}_R1_fastqc.zip */fastqc_data.txt > \
    ${qc_stats_dir}/${copa}_R1.pretrim.fastqc.data.txt
  unzip -p ${pretrim_qc_dir}/${copa}_R2_fastqc.zip */fastqc_data.txt > \
    ${qc_stats_dir}/${copa}_R2.pretrim.fastqc.data.txt
  unzip -p ${posttrim_qc_dir}/${copa}_R1.trim_fastqc.zip */fastqc_data.txt > \
    ${qc_stats_dir}/${copa}_R1.posttrim.fastqc.data.txt
  unzip -p ${posttrim_qc_dir}/${copa}_R2.trim_fastqc.zip */fastqc_data.txt > \
    ${qc_stats_dir}/${copa}_R2.posttrim.fastqc.data.txt
done < $file_index

echo "FastQC Results Unzipped"

seq_qual_base=${finalqc_dir}/sequence.quality.per.base.txt
echo -e "Sample\tPhase\tRead\tBase\tMean\tMedian\tLower Quartile\tUpper Quartile\t10th Percentile\t90th Percentile" > \
  $seq_qual_base
while IFS=, read $file_index_header
do
  # Skip the header row and any samples that should be excluded:
  [ "${Uid}" == "UID" ] && continue
  echo "Processing sample ${copa}"

  cat ${qc_stats_dir}/${copa}_R1.pretrim.fastqc.data.txt | \
    awk -v sample=$copa 'BEGIN { OFS="\t" }
      /^>>Per base sequence quality/ { flag=1; next }
      /^>>END_MODULE/ { flag=0 }
      (flag && $0 !~ /^#/) { print sample,"Pre-trim","R1",$0 }' >> \
      $seq_qual_base
  cat ${qc_stats_dir}/${copa}_R2.pretrim.fastqc.data.txt | \
    awk -v sample=$copa 'BEGIN { OFS="\t" }
      /^>>Per base sequence quality/ { flag=1; next }
      /^>>END_MODULE/ { flag=0 }
      (flag && $0 !~ /^#/) { print sample,"Pre-trim","R2",$0 }' >> \
      $seq_qual_base
  cat ${qc_stats_dir}/${copa}_R1.posttrim.fastqc.data.txt | \
    awk -v sample=$copa 'BEGIN { OFS="\t" }
      /^>>Per base sequence quality/ { flag=1; next }
      /^>>END_MODULE/ { flag=0 }
      (flag && $0 !~ /^#/) { print sample,"Post-trim","R1",$0 }' >> \
      $seq_qual_base
  cat ${qc_stats_dir}/${copa}_R2.posttrim.fastqc.data.txt | \
    awk -v sample=$copa 'BEGIN { OFS="\t" }
      /^>>Per base sequence quality/ { flag=1; next }
      /^>>END_MODULE/ { flag=0 }
      (flag && $0 !~ /^#/) { print sample,"Post-trim","R2",$0 }' >> \
      $seq_qual_base
done < $file_index

echo "Sequence per base extracted"

# Sequence duplication levels:
seq_dup_level=${finalqc_dir}/sequence.duplication.level.txt
echo -e "Sample\tPhase\tRead\tDuplication Level\tPercentage of deduplicated\tPercentage of total" > \
  $seq_dup_level
while IFS=, read $file_index_header
do
  # Skip the header row and any samples that should be excluded:
  [ "${Uid}" == "UID" ] && continue
  echo "Processing sample ${copa}"

  cat ${qc_stats_dir}/${copa}_R1.pretrim.fastqc.data.txt | \
    awk -v sample=$copa 'BEGIN { OFS="\t" }
      /^>>Sequence Duplication Levels/ { flag=1; next }
      /^>>END_MODULE/ { flag=0 }
      (flag && $0 !~ /^#/) { print sample,"Pre-trim","R1",$0 }' >> \
      $seq_dup_level
  cat ${qc_stats_dir}/${copa}_R2.pretrim.fastqc.data.txt | \
    awk -v sample=$copa 'BEGIN { OFS="\t" }
      /^>>Sequence Duplication Levels/ { flag=1; next }
      /^>>END_MODULE/ { flag=0 }
      (flag && $0 !~ /^#/) { print sample,"Pre-trim","R2",$0 }' >> \
      $seq_dup_level
  cat ${qc_stats_dir}/${copa}_R1.posttrim.fastqc.data.txt | \
    awk -v sample=$copa 'BEGIN { OFS="\t" }
      /^>>Sequence Duplication Levels/ { flag=1; next }
      /^>>END_MODULE/ { flag=0 }
      (flag && $0 !~ /^#/) { print sample,"Post-trim","R1",$0 }' >> \
      $seq_dup_level
  cat ${qc_stats_dir}/${copa}_R2.posttrim.fastqc.data.txt | \
    awk -v sample=$copa 'BEGIN { OFS="\t" }
      /^>>Sequence Duplication Levels/ { flag=1; next }
      /^>>END_MODULE/ { flag=0 }
      (flag && $0 !~ /^#/) { print sample,"Post-trim","R2",$0 }' >> \
      $seq_dup_level
done < $file_index

echo "Sequence duplication levels extracted"

# GC content levels:
gc_content_level=${finalqc_dir}/gc.content.level.txt
echo -e "Sample\tPhase\tRead\tGC Content\tCount" > \
  $gc_content_level
while IFS=, read $file_index_header
do
  # Skip the header row and any samples that should be excluded:
  [ "${Uid}" == "UID" ] && continue
  echo "Processing sample ${copa}"

  cat ${qc_stats_dir}/${copa}_R1.pretrim.fastqc.data.txt | \
    awk -v sample=$copa 'BEGIN { OFS="\t" }
      /^>>Per sequence GC content/ { flag=1; next }
      /^>>END_MODULE/ { flag=0 }
      (flag && $0 !~ /^#/) { print sample,"Pre-trim","R1",$0 }' >> \
      $gc_content_level
  cat ${qc_stats_dir}/${copa}_R2.pretrim.fastqc.data.txt | \
    awk -v sample=$copa 'BEGIN { OFS="\t" }
      /^>>Per sequence GC content/ { flag=1; next }
      /^>>END_MODULE/ { flag=0 }
      (flag && $0 !~ /^#/) { print sample,"Pre-trim","R2",$0 }' >> \
      $gc_content_level
  cat ${qc_stats_dir}/${copa}_R1.posttrim.fastqc.data.txt | \
    awk -v sample=$copa 'BEGIN { OFS="\t" }
      /^>>Per sequence GC content/ { flag=1; next }
      /^>>END_MODULE/ { flag=0 }
      (flag && $0 !~ /^#/) { print sample,"Post-trim","R1",$0 }' >> \
      $gc_content_level
  cat ${qc_stats_dir}/${copa}_R2.posttrim.fastqc.data.txt | \
    awk -v sample=$copa 'BEGIN { OFS="\t" }
      /^>>Per sequence GC content/ { flag=1; next }
      /^>>END_MODULE/ { flag=0 }
      (flag && $0 !~ /^#/) { print sample,"Post-trim","R2",$0 }' >> \
      $gc_content_level
done < $file_index

echo "GC Content extracted"

################################################################################
######################   Section 3: Alignment statistics   #####################
################################################################################

# In this section, we aggregate the alignment statistics provided by Bowtie2
# into a form that is easily imported into R.

alignment_stats=${finalqc_dir}/bowtie2.alignment.stats.txt
echo -e "Sample\tReads\tPairConcordantUnique\tPairConcordantMulti\tPairConcordantUnaligned\tPairDiscordantUnique\tMateUnique\tMateMulti\tMateUnaligned\tOverallAlignmentRate" > \
  $alignment_stats
while IFS=, read $file_index_header
do
  # Skip the header row and any samples that should be excluded:
  [ "${Uid}" == "UID" ] && continue

  echo "Processing sample ${copa}"

  cat ${logs_dir}/alignment/${Uid}.alignment.err | awk -v sample=$copa \
    'BEGIN{ OFS="\t";
            reads=0;
            paired=0;
            pair_concordant_unique=0
            pair_concordant_multi=0
            pair_concordant_unaligned=0
            pair_discordant_unique=0

            mate_unique=0
            mate_multi=0
            mate_unaligned=0
          }
     /^[[:digit:]]+ reads; of these:$/ { reads=$1 }
     /[[:digit:]]+ \([[:digit:]]+\.[[:digit:]]+\%\) were paired; of these:$/ { paired=$1 }
     /[[:digit:]]+ \([[:digit:]]+\.[[:digit:]]+\%\) aligned concordantly exactly 1 time$/ { pair_concordant_unique=$1 }
     /[[:digit:]]+ \([[:digit:]]+\.[[:digit:]]+\%\) aligned concordantly >1 times$/ { pair_concordant_multi=$1 }
     /[[:digit:]]+ \([[:digit:]]+\.[[:digit:]]+\%\) aligned concordantly 0 times$/ { pair_concordant_unaligned=$1 }
     /[[:digit:]]+ \([[:digit:]]+\.[[:digit:]]+\%\) aligned discordantly 1 time$/ { pair_discordant_unique=$1 }
     /[[:digit:]]+ \([[:digit:]]+\.[[:digit:]]+\%\) aligned exactly 1 time$/ { mate_unique=$1 }
     /[[:digit:]]+ \([[:digit:]]+\.[[:digit:]]+\%\) aligned >1 times$/ { mate_multi=$1 }
     /[[:digit:]]+ \([[:digit:]]+\.[[:digit:]]+\%\) aligned 0 times$/ { mate_unaligned=$1 }

     END{ total=pair_concordant_unique+pair_concordant_multi+pair_concordant_unaligned
          alignment_rate=(pair_concordant_unique+pair_concordant_multi+pair_discordant_unique+(mate_unique+mate_multi)/2)/total
          print sample,reads,pair_concordant_unique,pair_concordant_multi,pair_concordant_unaligned,pair_discordant_unique,mate_unique,mate_multi,mate_unaligned,alignment_rate
        }' >> \
    $alignment_stats
done < $file_index

echo "Alignment performance extracted" 

################################################################################
###################   Section 4: Post-alignment statistics   ###################
################################################################################

# In this section, we aggregate alignment summary metrics from various stages of
# filtering into a single file for subsequent import into R.

# Create an associative array containing the file suffixes for each stage:
declare -A stat_file_suffixes
stat_file_suffixes["Alignment"]="alignment.stats.txt"
stat_file_suffixes["Final"]="final.alignment.stats.txt"

# Collect picard alignment statistics for each stage of filtering:
picard_alignment=${finalqc_dir}/picard.alignment.stats.txt
echo -e "Sample\tPhase\tClass\tCount" > \
  $picard_alignment
while IFS=, read $file_index_header
do
  # Skip the header row and any samples that should be excluded:
  [ "$Uid" == "UID" ] && continue
  echo "Processing sample ${copa}"

  for stage in "${!stat_file_suffixes[@]}"; do
    cat ${alignment_stats_dir}/${copa}.${stat_file_suffixes[$stage]} | \
    awk -v sample=${copa} -v stage=${stage} \
      'BEGIN{ OFS="\t" }
       /^PAIR/ { print sample,stage,"Total Reads",$2
                 print sample,stage,"PF Reads",$3
                 print sample,stage,"PF Reads Aligned",$6
                 print sample,stage,"PF HQ Aligned Reads",$9
                 print sample,stage,"Reads Aligned in Pairs",$17
               }' >> $picard_alignment
  done
done < $file_index

unset stat_file_suffixes

echo "Picard evaluation extracted"

    # alignment stats
    # last row is paired QC
    # total reads, PF, reads, PF reads aligned, PF HQ aligned bases, PF HQ aligned bases Q>20, ave read length, paired reads aligned
#    preshiftstats=$( sed -n '7p;10p' ${alignment_stats_dir}/${uid}.alignment.stats.txt | cut -f2,3,6,10,11,16,17 )

    # alignment stats after shifting reads
    # last row is paired QC
    # total reads, PF, reads, PF reads aligned, PF HQ aligned bases, PF HQ aligned bases Q>20, ave read length, paired reads aligned
 #   postshiftstats=$( sed -n '7p;10p' ${alignment_stats_dir}/${uid}.alignment.stats.shifted.txt | cut -f2,3,6,10,11,16,17 )

    # duplication rates, histogram shows how much greater coverage to expect from x times many actual coverage. Should plateau at a point.
#     dupstats=$( sed -n '7,8p' ${logs_dir}/${uid}.atac.duplicates.log | cut --complement -f1 )

#     # number of peaks, ave peak width
#     peaknum=$( wc -l ${peaks_dir}/${uid}_peaks.narrowPeak | cut -f1 -d" " )
#     ave_width=$( Rscript ${src_dir}/ave_peakwidth.R ${peaks_dir}/${uid}_peaks.narrowPeak | cut -f2 -d" " )

#     # Num of fp discovered, num of peaks w/fp, num of unique TFs with fp, num 
#     # num fp
#     fps=$( wc -l ${fp_dir}/${uid}.bed | cut -f1 -d" " )
#     # num peaks w/fp
#     pfps=$( intersectBed -u -a ${peaks_dir}/${uid}_peaks.narrowPeak -b ${fp_dir}/${uid}.bed | wc -l | cut -f1 -d" " )
#     # num unique human pwms w/fps found
#     unique_pwms=$( cut -f4 ${match_dir}/${uid}_mpbs.bed | sort | grep "HUMAN" | sort | uniq | wc -l )
#     # num unique human tfs w/fps found
#     unique_tfs=$( cut -f4 ${match_dir}/${uid}_mpbs.bed | sort | grep "HUMAN" | cut -f1 -d_ | uniq | wc -l )

#     echo "Preshift: ${preshiftstats}
# Postshift: ${postshiftstats} 
# Duplication Stats: ${dupstats} 
# Number of peaks: ${peaknum} 
# Average peakwidth: ${ave_width} 
# Footprints discovered: ${fps} 
# Peaks with footprints: ${pfps} 
# Unique Human PWMs found: ${unique_pwms}
# Unique Human TFs found: ${unique_tfs}" > ${results_dir}/${uid}_metrics.txt
