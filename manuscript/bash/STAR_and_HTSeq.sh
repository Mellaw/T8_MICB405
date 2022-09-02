#!/usr/bin/env bash

# MICB 405: Group 8
# STAR (RNA alignment) and HTSeq (determining reads/gene)
# Dec 7, 2021

mkdir /projects/micb405/analysis/GROUP8/star && cd /projects/micb405/analysis/GROUP8/star
mkdir STARIndex

### STAR ###

#generate STAR index using the mm10 reference genome (.fasta) and gene file (.gtf)
STAR \
  --runMode genomeGenerate \
  --genomeDir STARIndex \
  --genomeFastaFiles /projects/micb405/analysis/STAR_tutorial/Mus_musculus.GRCm38.dna.primary_assembly.fa \
  --sjdbGTFfile /projects/micb405/analysis/STAR_tutorial/Mus_musculus.GRCm38.84.gtf \
  --sjdbOverhang 49 \
  --runThreadN 16

#align RNA-seq data (.fastq) to your STAR index
for FILE in /projects/micb405/analysis/GROUP8/fastq/*
do
 NAME=$(basename $FILE)
 STAR \
  --genomeDir STARIndex/ \
  --readFilesIn $FILE \
  --outFileNamePrefix /projects/micb405/analysis/GROUP8/RNA_alignments2/$NAME \
  --runThreadN 8 \
  --limitBAMsortRAM 60000000000 \
  --outSAMattrRGline ID:$NAME SM:$NAME \
  --outBAMsortingThreadN 8 \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMunmapped Within \
  --outSAMstrandField intronMotif \
  --readFilesCommand zcat \
  --chimSegmentMin 20 \
  --genomeLoad NoSharedMemory
done

### HTSeq ###

mkdir /projects/micb405/analysis/GROUP8/htseq_output && cd /projects/micb405/analysis/GROUP8/htseq_output

#install HTSeq
pip3 install htseq

#index your STAR output .bam files (HTSeq requirement)
#keep indexed .bai files in the same directory as the original .bam files
samtools index /projects/micb405/analysis/GROUP8/RNA_alignments2/*.bam

#run HTSeq on your STAR output .bam files to determine reads/gene
for FILE in /projects/micb405/analysis/GROUP8/RNA_alignments2/*.fastq.gzAligned.sortedByCoord.out.bam
do
 NAME=$(basename $FILE)
htseq-count \
  $FILE \
  /projects/micb405/analysis/STAR_tutorial/Mus_musculus.GRCm38.84.gtf \
  -f bam \
  -r pos \
  --stranded=no \
  > ${NAME}.htseq.out
done

#export HTSeq files
scp hsokolovska_mb21@orca2.bcgsc.ca:/projects/micb405/analysis/GROUP8/htseq_output/* /Users/Helena/405_Linux_outputs
