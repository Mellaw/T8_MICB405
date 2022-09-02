# T8_MICB405

Sunset Enhancers: Tracing H3K27 Acetylation on Closed Chromatin in Myeloid Lineage Differentiation

This is the full code repository for Group 8's MICB 405 capstone project. Python and R were used for data manipulation and figure generation. All iChIP-seq and ATAC-seq pre-processing is done in bash. 

After downloading this repository, follow these step to replicate our findings:

To download, process, and analyze ATAC- and iCHIP-seq data, run ./manuscript/bash/1_Data_Import_and_ChIP-seq.sh **from the ./manuscript folder**. Figures can be recreated using the data generated (final data located in data/enhancer_status.tsv), but exact code is not included. The following dependencies are required:
- Ubuntu 16.04.5
- bwa 0.7.17-r1188
- bigBedToBed v1
- fastq-dump 2.9.2
- GNU Parallel
- samtools 1.9 (using htslib 1.9)
- macs2 2.2.7.1
- bedtools v2.27.1
- Python 3.9.7 (with the follow packages): These 3 packages will take care of the dependencies
  - pandas 1.3.3
  - numpy 1.20.3
  - seaborn 0.11.2

For the full pipeline of analysis, run numbered scripts in order:

1_Data_Import_and_ChIP-seq.sh
   Calls on scripts in python folder:
     get_desired_runs.py
     rename_fastq.py
     get_median_enh_treat_control_pileup.py
     get_max_enhancer_q_value.py
2_RNA-seq_STAR_and_HTSeq.sh
3_Assign_Enhancer_Status.py
4_Calculate_TPM.py
5_DESeq2_and_GO_analysis

Contributing Authors: Andy Murtha, Helena Sokolovska, Kitty Martens, Melanie Law, Kalen Dofher, and Annice Li
