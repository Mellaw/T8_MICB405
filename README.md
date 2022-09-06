# T8_MICB405

Sunset Enhancers: Tracing H3K27 Acetylation on Closed Chromatin in Myeloid Lineage Differentiation

Contributing Authors: Melanie Law, Helena Sokolovska, Andy Murtha, Kitoosepe Martens, Annice Li, and Kalen Dofher

This is the full code repository for our study. Python and R were used for data manipulation and figure generation. All iChIP-seq and ATAC-seq pre-processing is done in bash. 

After downloading this repository, follow these step to replicate our findings:

To download, process, and analyze ATAC- and iCHIP-seq data, run ./manuscript/bash/1_Data_Import_and_ChIP-seq.sh **from the ./manuscript folder**. Figures can be recreated using the data generated (final data located in data/enhancer_status.tsv), but exact code is not included.

The following dependencies are required:
- Ubuntu 16.04.5
- bwa 0.7.17-r1188
- bigBedToBed v1
- fastq-dump 2.9.2
- GNU Parallel
- samtools 1.9 (using htslib 1.9)
- macs2 2.2.7.1
- bedtools v2.27.1
- Python 3.9.7 (with the following packages): These 3 packages will take care of the dependencies
  - pandas 1.3.3
  - numpy 1.20.3
  - seaborn 0.11.2
- Genomic Regions Enrichment of Annotations Tool (GREAT) 4.0.4
- STAR 2.7.9a
- HTSeq 0.11.2
- RStudio 2021.09.0
  - DESeq2 1.30.1 
  - dplyr 1.0.8
  - tidyverse 1.3.1 
  - xlsx 0.6.5
  - AnnotationDbi 1.52.0
  - org.Mm.eg.db 3.12.0  
  - GO.db 3.12.1 
  - GOstats 2.56.0

For the full pipeline of analysis, run numbered scripts in order:
- 1_Data_Import_and_ChIP-seq.sh
  - Calls on scripts in python folder:
     - get_desired_runs.py
     - rename_fastq.py
     - get_median_enh_treat_control_pileup.py
     - get_max_enhancer_q_value.py
- 2_RNA-seq_STAR_and_HTSeq.sh
- 3_Assign_Enhancer_Status.py
- 4_Calculate_TPM.py
- 5_DESeq2_and_GO_analysis
