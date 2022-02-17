#!/bin/bash

REFERENCE=ZIKV-PR.fasta
GFF_FILE=ZIKV-PR.gff3
PRIMERS=primers.fasta
PRIMER_PAIRS=primer_pair_information.tsv
THRESHOLD=0.03

date
echo "Mutation detection threshold, percent:"
awk "BEGIN {print $THRESHOLD*100}"
# We are counting original BAM files here because one extra BAM file will be created soon for primers.
TOTAL_FILES=`find -maxdepth 1 -iname '*.bam' | wc -l`
ARR=($(ls *.bam))
bwa index $REFERENCE
bwa mem -k 5 -T 16 $REFERENCE $PRIMERS | samtools view -b -F 4 > primers.bam
bedtools bamtobed -i primers.bam > primers.bed
printf "\n\n"     
echo "Part 1: processing individual files with indexing, trimming, sorting, and indexing"
printf "\n"
for ((i=0; i<$TOTAL_FILES; i+=1)) {
    FILE_NAME=`echo ${ARR[$i]} | awk -F "." '{print $1}'`
    printf "\n"
    echo "[processing] $FILE_NAME"
    echo "[indexing] samtools index $FILE_NAME.bam"
    samtools index $FILE_NAME.bam
    echo "[trimming] ivar trim -b primers.bed -i ${ARR[$i]} -p $FILE_NAME.trimmed"
    ivar trim -b primers.bed -i ${ARR[$i]} -p $FILE_NAME.trimmed
    echo "[sorting] samtools sort $FILE_NAME.trimmed.bam -o $FILE_NAME.trimmed.sorted.bam"
    samtools sort $FILE_NAME.trimmed.bam -o $FILE_NAME.trimmed.sorted.bam
    echo "[indexing] samtools index $FILE_NAME.trimmed.sorted.bam"
    samtools index $FILE_NAME.trimmed.sorted.bam
}
printf "\n\n"
echo "Part 2: processing two replicates with merging, making and indexing consensus, primers processing (aligning, framing, and finding mismatches)"
printf "\n"
for ((i=1; i<$TOTAL_FILES; i+=2)) {
    SAMPLE_NAME=`echo ${ARR[$i]} | awk -F "_" '{print $1}'`
    FILE_NAME1=`echo ${ARR[$i-1]} | awk -F "." '{print $1}'`
    FILE_NAME2=`echo ${ARR[$i]} | awk -F "." '{print $1}'`
    printf "\n"
    echo "[processing] $SAMPLE_NAME"
    echo "[merging replicates] samtools merge $SAMPLE_NAME.merged.bam $FILE_NAME1.trimmed.sorted.bam $FILE_NAME2.trimmed.sorted.bam"
    samtools merge $SAMPLE_NAME.merged.bam $FILE_NAME1.trimmed.sorted.bam $FILE_NAME2.trimmed.sorted.bam
    echo "[making consensus] samtools mpileup -A -d 1000000 -Q 0 $SAMPLE_NAME.merged.bam | ivar consensus -p $SAMPLE_NAME.consensus"
    samtools mpileup -A -d 1000000 -Q 0 $SAMPLE_NAME.merged.bam | ivar consensus -p $SAMPLE_NAME.consensus
    echo "[primers processing: indexing consensus] bwa index $SAMPLE_NAME.consensus.fa"
    bwa index $SAMPLE_NAME.consensus.fa
    echo "[primers processing: conversion to sorted bam] bwa mem -k 5 -T 16 $SAMPLE_NAME.consensus.fa $PRIMERS | samtools view -bS -F 4 | samtools sort -o $SAMPLE_NAME.primers_consensus.bam"
    bwa mem -k 5 -T 16 $SAMPLE_NAME.consensus.fa $PRIMERS | samtools view -bS -F 4 | samtools sort -o $SAMPLE_NAME.primers_consensus.bam
    echo "[primers processing: creating BED] bedtools bamtobed -i $SAMPLE_NAME.primers_consensus.bam > $SAMPLE_NAME.primers_consensus.bed"
    bedtools bamtobed -i $SAMPLE_NAME.primers_consensus.bam > $SAMPLE_NAME.primers_consensus.bed
    echo "[primers processing: creating TSV file] samtools mpileup -A -d 1000000 --reference $SAMPLE_NAME.consensus.fa -Q 0 $SAMPLE_NAME.primers_consensus.bam | ivar variants -p $SAMPLE_NAME.primers_consensus -t $THRESHOLD"
    samtools mpileup -A -d 1000000 --reference $SAMPLE_NAME.consensus.fa -Q 0 $SAMPLE_NAME.primers_consensus.bam | ivar variants -p $SAMPLE_NAME.primers_consensus -t $THRESHOLD
    echo "[primers processing: masking] ivar getmasked -i $SAMPLE_NAME.primers_consensus.tsv -b $SAMPLE_NAME.primers_consensus.bed -f $PRIMER_PAIRS -p $SAMPLE_NAME.mismatches"
    ivar getmasked -i $SAMPLE_NAME.primers_consensus.tsv -b $SAMPLE_NAME.primers_consensus.bed -f $PRIMER_PAIRS -p $SAMPLE_NAME.mismatches
}
printf "\n\n"
echo "Part 3: processing individual files with removing reads with mismatches, sorting BAM-file, and calling SNVs"
printf "\n"
for ((i=0; i<$TOTAL_FILES; i+=1)) {
    SAMPLE_NAME=`echo ${ARR[$i]} | awk -F "_" '{print $1}'`
    FILE_NAME=`echo ${ARR[$i]} | awk -F "." '{print $1}'`
    printf "\n"
    echo "[processing] $FILE_NAME"
    echo "[removing reads] ivar removereads -i $FILE_NAME.trimmed.sorted.bam -p $FILE_NAME.trimmed.sorted.masked.bam -t $SAMPLE_NAME.mismatches.txt -b primers.bed"
    ivar removereads -i $FILE_NAME.trimmed.sorted.bam -p $FILE_NAME.trimmed.sorted.masked.bam -t $SAMPLE_NAME.mismatches.txt -b primers.bed
    echo "[sorting resulting BAM] samtools sort $FILE_NAME.trimmed.sorted.masked.bam -o $FILE_NAME.trimmed.sorted.masked.sorted.bam"
    samtools sort $FILE_NAME.trimmed.sorted.masked.bam -o $FILE_NAME.trimmed.sorted.masked.sorted.bam
    echo "[generating final TSV file] samtools mpileup -A -d 1000000 --reference $REFERENCE -B -Q 0 $FILE_NAME.trimmed.sorted.masked.sorted.bam | ivar variants -p $FILE_NAME.SNV -t $THRESHOLD -r $REFERENCE -g GFF_FILE"
    samtools mpileup -A -d 1000000 --reference $REFERENCE -B -Q 0 $FILE_NAME.trimmed.sorted.masked.sorted.bam | ivar variants -p $FILE_NAME.SNV -t $THRESHOLD -r $REFERENCE -g $GFF_FILE
}
printf "\n\n"
echo "Part 4: Filtering the same SNVs from replicates"
printf "\n"
for ((i=1; i<$TOTAL_FILES; i+=2)) {
    SAMPLE_NAME=`echo ${ARR[$i]} | awk -F "_" '{print $1}'`
    FILE_NAME1=`echo ${ARR[$i-1]} | awk -F "." '{print $1}'`
    FILE_NAME2=`echo ${ARR[$i]} | awk -F "." '{print $1}'`
    printf "\n"
    echo "[processing] $SAMPLE_NAME"
    echo "[filtering SNVs across replicates] ivar filtervariants -p $SAMPLE_NAME.SNV_final $FILE_NAME1.SNV.tsv $FILE_NAME2.SNV.tsv"
    ivar filtervariants -p $SAMPLE_NAME.SNV_final $FILE_NAME1.SNV.tsv $FILE_NAME2.SNV.tsv
}
echo "[Sorting files to directories and renaming]"
mkdir 0.input && mv $PRIMERS $REFERENCE $PRIMER_PAIRS $GFF_FILE $_
mkdir 1.index_files && mv *.pac $_ && mv *.sa $_ && mv *.ann $_ && mv *.amb $_ && mv *.fai $_ && mv *.bai $_ && mv *.bwt $_
mkdir 1.Primers_info && mv primers.bam primers.bed $_
mkdir 2a.Consensus && mv *.consensus.fa $_
mkdir 2c.BED && mv *.bed $_
mkdir 2d.Consensus-TSV && mv *.primers_consensus.tsv $_
mkdir 3.Mismatches  && mv *.mismatches.txt $_
mkdir 4.Quality && mv *.qual.txt $_
mkdir 7.SNV-final && mv *.SNV_final.tsv $_
mkdir 6.SNV-singlets && mv *.SNV.tsv $_
mkdir 5.BAM-final && mv *.masked.sorted.bam $_
rm *.masked.bam
mkdir 2b.BAM-primers && mv *.primers_consensus.bam $_
rm *.merged.bam
rm *.sorted.bam
rm *.trimmed.bam
mkdir 0.BAM-original && mv *.bam $_
cd ./7.SNV-final
for FILE in *; do mv "$FILE" "${FILE/SNV_final./}"; done
cd ../5.BAM-final
for FILE in *; do mv "$FILE" "${FILE/trimmed.sorted.masked.sorted./}"; done
cd ..
date
echo "[FINISHED]"
