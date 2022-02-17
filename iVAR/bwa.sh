#!/bin/bash

INDEX=illumina-phix.fasta
#INDEX=BaPCV2b.fasta
#INDEX=ZIKV-PR.fasta
SEARCH_PATTERN='*.fastq.gz'

date
echo "[building index] $INDEX"
bwa index $INDEX
TOTAL_FILES=`find ./ -maxdepth 1 -iname "*$SEARCH_PATTERN*" | wc -l`
ARR=($(ls $SEARCH_PATTERN))
for ((i=0; i<$TOTAL_FILES; i+=2))
{
    printf "\n"
    SAMPLE_NAME=`echo ${ARR[$i]} | awk -F "_" '{print $1}'`
    echo "[mapping running for $[i/2+1]/$[TOTAL_FILES/2]] $SAMPLE_NAME"
    echo "bwa mem -t 32 $INDEX ${ARR[$i]} ${ARR[$i+1]} | samtools view -b -F 4 -F 2048 | samtools sort -o $SAMPLE_NAME.bam"
    bwa mem -t 32 $INDEX ${ARR[$i]} ${ARR[$i+1]} | samtools view -b -F 4 -F 2048 | samtools sort -o $SAMPLE_NAME.bam
}
date