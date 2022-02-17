#!/bin/bash

mkdir 1-merged
mkdir 2-sorted
TOTAL_FILES=`find ./ -maxdepth 1 -iname '*.bam' | wc -l`
ARR=($(ls *.bam))
for ((i=0; i<$TOTAL_FILES; i+=2))
{
    SAMPLE_NAME=`echo ${ARR[$i]} | awk -F "_" '{print $1}'`
    printf "\n"
    echo "[merging] $SAMPLE_NAME"
    echo "samtools merge ./1-merged/$SAMPLE_NAME.bam ${ARR[$i]} ${ARR[$i+1]}"
    samtools merge ./1-merged/$SAMPLE_NAME.bam ${ARR[$i]} ${ARR[$i+1]}
}
cd 1-merged
TOTAL_FILES=`find ./ -maxdepth 1 -iname '*.bam' | wc -l`
ARR=($(ls *.bam))
for ((i=0; i<$TOTAL_FILES; i+=1))
{
    SAMPLE_NAME=`echo ${ARR[$i]} | awk -F "_" '{print $1}'`
    printf "\n"
    echo "[sorting] $SAMPLE_NAME "
    echo "samtools view -b ${ARR[$i]} | samtools sort -o ./../2-sorted/$SAMPLE_NAME"
    samtools view -b ${ARR[$i]} | samtools sort -o ./../2-sorted/$SAMPLE_NAME
}
cd ..
