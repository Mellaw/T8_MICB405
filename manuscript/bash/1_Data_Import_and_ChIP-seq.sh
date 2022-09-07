#!/usr/bin/env bash

#Sunset Enhancers: Tracing H3K27 Acetylation on Closed Chromatin in Myeloid Lineage Differentiation
#Authors: Melanie Law, Helena Sokolovska, Andrew Murtha, Kitoosepe Martens, Annice Li, and Kalen Dofher
#Date: Sept 5, 2022

##### DATA IMPORT ####

wget https://hgdownload.soe.ucsc.edu/gbdb/mm10/encode3/encode3RenInteractEnhancerAll.bb -P ./data/references
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz -P ./data/references

gunzip ./data/references/mm10.fa.gz
bigBedToBed ./data/references/encode3RenInteractEnhancerAll.bb ./data/references/mm10_enhancerAll.bed

python ./srr_metadata/get_desired_runs.py

metadata=./srr_metadata/srr_to_download.tsv
outdir=./data/fastq/

rm ./srr_metadata/gnu_download.txt
touch ./srr_metadata/gnu_download.txt

i=1
for srr in $(tail -n +2 $metadata | awk '{print $1}');
do
        echo $srr;
        echo "fastq-dump --gzip --outdir $outdir $srr" >> ./srr_metadata/gnu_download.txt
done;

parallel -j16 < ./srr_metadata/gnu_download.txt

#### ALIGN BAM FILES FOR iCHIP AND ATAC-SEQ #####

python ./python/rename_fastq.py


fq_dir=./data/fastq
ref=/projects/micb405/analysis/references/mm10/mm10.fa
outdir=./data/DNA_alignments

pfile=./parallel/alignment.txt
rm $pfile
touch $pfile

i=1
for f in ./data/fastq/ATAC-Seq_*.fastq.gz;
do
        echo $f;
        echo "bwa mem -t2 $ref $fq_dir/$f | samtools view -b -o $outdir/${f/.fastq.gz/.bam};" >> $pfile
done;

for f in H3K27Ac_*.fastq.gz;
do
        echo $f;
        echo "bwa mem -t2 $ref $fq_dir/$f | samtools view -b -o $outdir/${f/.fastq.gz/.bam};" >> $pfile
done;

for f in H3K4me1_*.fastq.gz;
do
        echo $f;
        echo "bwa mem -t2 $ref $fq_dir/$f | samtools view -b -o $outdir/${f/.fastq.gz/.bam};" >> $pfile
done;

parallel -j12 < $pfile

#### RUN MACS2 ON ALL iCHIP AND ATAC-SEQ DATA #####


# cd into alignments folder
bam_dir=./data/DNA_alignments

mkdir ./data/macs2_comb/
touch ./parallel/macs2_parallel.txt

# Loop over all bam files and write a macs2 call to the parallel file
for f in $(echo $bam_dir/*.bam | tr " " "\n" | cut -d"_" -f1,2 | uniq);
do
	fs=$(echo $bam_dir/${f}*)
    echo "macs2 callpeak -t ${fs} -g mm -n ${f/.bam/} -B --outdir ./data/macs2_comb/" >> ./parallel/macs2_parallel.txt;
done
parallel -j12 < ./parallel/macs2_parallel.txt

#### CREATE FOLD-ENRICHMENT FILES FROM MACS2 BEDGRAPH OUTPUT ####

#cd into macs2 folder
macs_dir=./data/macs2_comb

rm ../parallel/get_FE.txt
touch ../parallel/get_FE.txt

#loop over types
for f in $macs_dir/*treat_pileup.bdg; do
	t=$(echo $f | cut -d"_" -f1,2)
	treat=$macs_dir/${f}
	ctrl=$macs_dir/${t}_control_lambda.bdg
	outfe=$macs_dir/${t}_FE.bdg
	outq=$macs_dir/${t}_q.bdg
	echo ${treat},${ctrl},${outfe}
	echo "macs2 bdgcmp -t $treat -c $ctrl -o $outfe -m FE" >> ../parallel/get_FE.txt
	echo "macs2 bdgcmp -t $treat -c $ctrl -o $outq -m qpois" >> ../parallel/get_FE.txt
done;

parallel -j12 < ../parallel/get_FE.txt

#### SORT FE BEDGRAPH FILES

cd $macs_dir/
parallel_dir=../../parallel
rm $parallel_dir/sort_FE.txt
touch $parallel_dir/sort_FE.txt

#loop over types
for f in *FE.bdg; do
        echo "bedtools sort -i $f > ${f/.bdg/.sorted.bdg}" >> $parallel_dir/sort_FE.txt
done;
for f in *q.bdg; do
        echo "bedtools sort -i $f > ${f/.bdg/.sorted.bdg}" >> $parallel_dir/sort_FE.txt
done;

parallel -j12 < $parallel_dir/sort_FE.txt

for f in *FE.sorted.bdg; do
    mv $f ${f/.sorted.bdg/.bdg}
done;

for f in *q.sorted.bdg; do
        mv $f ${f/.sorted.bdg/.bdg}
done;


### Merge biological duplicates ####

enh=../references/mm10_enhancerAll.bed

rm $parallel_dir/merge_enh.txt
touch $parallel_dir/merge_enh.txt

#loop over types
for f in *_FE.bdg; do
        t=$(echo $f | cut -d"_" -f1,2);
        echo $t
        a=$(echo $t*FE.bdg | cut -d" " -f1);
        echo "bedtools intersect -wa -wb -a $enh -b $a -filenames -sorted > ../comb_bdg/${t}_FE.bdg"
        echo "bedtools intersect -wa -wb -a $enh -b $a -filenames -sorted > ../comb_bdg/${t}_FE.bdg" >> $parallel_dir/merge_dups.txt
done

for f in *_q.bdg; do
        t=$(echo $f | cut -d"_" -f1,2);
        echo $t
        echo "bedtools intersect -wa -wb -a $enh -b $f -filenames -sorted > ../comb_bdg/${t}_q.bdg"
        echo "bedtools intersect -wa -wb -a $enh -b $f -filenames -sorted > ../comb_bdg/${t}_q.bdg" >> $parallel_dir/merge_dups.txt
done

parallel -j12 < ../parallel/merge_dups.txt

## Merge peak files

rm ../parallel/merge_peaks.txt
touch ../parallel/merge_peaks.txt

#loop over types
for f in *_peaks.narrowPeak; do
	t=$(echo $f | cut -d"_" -f1,2);
	echo $t
	echo "bedtools intersect -wa -wb -a $enh -b $f -filenames -sorted > ../comb_bdg/${t}_peaks.bdg"
	echo "bedtools intersect -wa -wb -a $enh -b $f -filenames -sorted > ../comb_bdg/${t}_peaks.bdg" >> ../parallel/merge_peaks.txt
done

parallel -j12 < ../parallel/merge_peaks.txt
cd ../../
python ./python/get_median_enh_treat_control_pileup.py
python ./python/get_max_enhancer_q_value.py
