#!/usr/bin/env bash

##### DATA IMPORT ####
rm /projects/micb405/analysis/GROUP8/fastq/*

python /projects/micb405/analysis/GROUP8/srr_metadata/get_desired_runs.py

metadata=/projects/micb405/analysis/GROUP8/srr_metadata/srr_to_download.tsv
outdir=/projects/micb405/analysis/GROUP8/fastq/

rm /projects/micb405/analysis/GROUP8/srr_metadata/gnu_download.txt
touch /projects/micb405/analysis/GROUP8/srr_metadata/gnu_download.txt

i=1
for srr in $(tail -n +2 $metadata | awk '{print $1}');
do
        echo $srr;
        echo "fastq-dump --gzip --outdir $outdir $srr" >> /projects/micb405/analysis/GROUP8/srr_metadata/gnu_download.txt
done;

parallel -j16 < /projects/micb405/analysis/GROUP8/srr_metadata/gnu_download.txt

#### ALIGN BAM FILES FOR iCHIP AND ATAC-SEQ #####

python /projects/micb405/analysis/GROUP8/srr_metadata/rename_fastq.py


cd /projects/micb405/analysis/GROUP8/fastq

ref=/projects/micb405/analysis/references/mm10/mm10.fa
outdir=/projects/micb405/analysis/GROUP8/DNA_alignments

pfile=/projects/micb405/analysis/GROUP8/parallel/alignment.txt
rm $pfile
touch $pfile

i=1
for f in ATAC-Seq_*.fastq.gz;
do
        echo $f;
        echo "bwa mem -t2 $ref $f | samtools view -b -o $outdir/${f/.fastq.gz/.bam};" >> $pfile
done;

for f in H3K27Ac_*.fastq.gz;
do
        echo $f;
        echo "bwa mem -t2 $ref $f | samtools view -b -o $outdir/${f/.fastq.gz/.bam};" >> $pfile
done;

for f in H3K4me1_*.fastq.gz;
do
        echo $f;
        echo "bwa mem -t2 $ref $f | samtools view -b -o $outdir/${f/.fastq.gz/.bam};" >> $pfile
done;

parallel -j12 < $pfile

#### RUN MACS2 ON ALL iCHIP AND ATAC-SEQ DATA #####


# cd into alignments folder
cd /projects/micb405/analysis/GROUP8/DNA_alignments

rm -r ../macs2_comb/*

rm ../parallel/macs2_parallel.txt
touch ../parallel/macs2_parallel.txt

# Loop over all bam files and write a macs2 call to the parallel file
for f in $(echo *.bam | tr " " "\n" | cut -d"_" -f1,2 | uniq);
do
	fs=$(echo ${f}*)
    echo "macs2 callpeak -t ${fs} -g mm -n ${f/.bam/} -B --outdir ../macs2_comb/" >> ../parallel/macs2_parallel.txt;
done
parallel -j12 < /projects/micb405/analysis/GROUP8/parallel/macs2_parallel.txt

#### CREATE FOLD-ENRICHMENT FILES FROM MACS2 BEDGRAPH OUTPUT ####

#cd into macs2 folder
cd /projects/micb405/analysis/GROUP8/macs2_comb

rm ../parallel/get_FE.txt
touch ../parallel/get_FE.txt

#loop over types
for f in *treat_pileup.bdg; do
	t=$(echo $f | cut -d"_" -f1,2)
	treat=${f}
	ctrl=${t}_control_lambda.bdg
	outfe=${t}_FE.bdg
	outq=${t}_q.bdg
	echo ${treat},${ctrl},${outfe}
	echo "macs2 bdgcmp -t $treat -c $ctrl -o $outfe -m FE" >> ../parallel/get_FE.txt
	echo "macs2 bdgcmp -t $treat -c $ctrl -o $outq -m qpois" >> ../parallel/get_FE.txt
done;

parallel -j12 < ../parallel/get_FE.txt

#### SORT FE BEDGRAPH FILES

cd /projects/micb405/analysis/GROUP8/macs2_comb

rm ../parallel/sort_FE.txt
touch ../parallel/sort_FE.txt

#loop over types
for f in *FE.bdg; do
        echo "bedtools sort -i $f > ${f/.bdg/.sorted.bdg}" >> ../parallel/sort_FE.txt
done;
for f in *q.bdg; do
        echo "bedtools sort -i $f > ${f/.bdg/.sorted.bdg}" >> ../parallel/sort_FE.txt
done;

cat ../parallel/sort_FE.txt

parallel -j12 < ../parallel/sort_FE.txt

for f in *FE.sorted.bdg; do
        mv $f ${f/.sorted.bdg/.bdg}
done;

for f in *q.sorted.bdg; do
        mv $f ${f/.sorted.bdg/.bdg}
done;


### Merge biological duplicates ####
cd /projects/micb405/analysis/GROUP8/macs2_comb

enh=/projects/micb405/analysis/GROUP8/mm10_enhancerAll.bed

rm ../parallel/merge_enh.txt
touch ../parallel/merge_enh.txt

#loop over types
for f in *_FE.bdg; do
        t=$(echo $f | cut -d"_" -f1,2);
        echo $t
        a=$(echo $t*FE.bdg | cut -d" " -f1);
        echo "bedtools intersect -wa -wb -a $enh -b $a -filenames -sorted > ../comb_bdg/${t}_FE.bdg"
        echo "bedtools intersect -wa -wb -a $enh -b $a -filenames -sorted > ../comb_bdg/${t}_FE.bdg" >> ../parallel/merge_dups.txt
done

for f in *_q.bdg; do
        t=$(echo $f | cut -d"_" -f1,2);
        echo $t
        echo "bedtools intersect -wa -wb -a $enh -b $f -filenames -sorted > ../comb_bdg/${t}_q.bdg"
        echo "bedtools intersect -wa -wb -a $enh -b $f -filenames -sorted > ../comb_bdg/${t}_q.bdg" >> ../parallel/merge_dups.txt
done

parallel -j12 < ../parallel/merge_dups.txt

## Merge peak files

cd /projects/micb405/analysis/GROUP8/macs2_comb

enh=/projects/micb405/analysis/GROUP8/mm10_enhancer.bed

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

python /projects/micb405/analysis/GROUP8/src/get_median_enh_treat_control_pileup.py
python /projects/micb405/analysis/GROUP8/src/get_max_enhancer_q_value.py




