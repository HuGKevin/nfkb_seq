This is a list of all the files in this directory and what they do. Bolded files are required for the pipelines to function.

* __atac_env_package_list.txt__ - List of packages in the miniconda environment "atac_env" used for read trimming. This can be used to duplicate that original environment. 
* __ATAC_MASTER.sh__ - Contains all the job array scripts necessary for running pipeline on ATAC-seq data. 
* __bind_clusters.R__ - Script used to convert MCL output to a bed file of cluster coordinates. 
* __dba_atac.R__ - Code for annotating clusters, running DiffBind to find differentially accessible clusters, and annotating those clusters. 
* downsample.sh - can be deleted, merged into master scripts. 
* fastqc_summary.sh - Script used to puill out from FastQC reports which tests each sequencing run passed, warninged, or failed, with the intention of using this to QC those individual runs. Ultimately didn't use, but will keep in case it's useful for something else. 
* figures.R - Scripts making figures for genomics analysis. For the most part, they're figures regarding the mcl clusters. 
* lib_depths.sh - Scripts extracting library depths from log files
* mcl_bias.R - Script creating plots to test whether certain libraries have more active singleton clusters than others. 
* __MINT_MASTER.sh__ - Contains all the job array scripts necessary for running pipeline on MintChIP data. 
* multiqc_analysis.R - Script analyzing output of multifastqc reports. This wasn't actually very useful though. 
* peak_cov.R - Script for generating table of what proportion of genome each library's peaks cover.
* __peak_matrix.R__ - Script for generating presence/absence matrix from MCL analysis output. Necessary for pipeline. Also contains some trash code from trying to hierarchically cluster libraries. 
* __pyadapter_trim.py__ - Script for trimming adapters from reads in fastq format. Necessary in pipeline.
* __shift_reads.py__ - Script for shifting aligned reads to correct for transposase bias. Necessary in ATAC-seq pipeline.
* test_seq_correlates.R - This was testing whether sequencing depth of libraries correlated with biological variables. I only tested whether depth was related to rs1800693 genotype using anovas and Kruskal-wallis tests, which seemed to give significant results. Further analysis warranted. 
