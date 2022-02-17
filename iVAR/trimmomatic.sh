#!/bin/bash
date
mkdir trimmed
FOLDER=/home/user/bin/Trimmomatic-0.39
TOTAL_FILES=`find -iname '*.gz' | wc -l`
ARR=($(ls *.gz))
for ((i=0; i<$TOTAL_FILES; i+=2))
{
    SAMPLE_NAME=`echo ${ARR[$i]} | awk -F "_" '{print $1}'`
    printf "\n"
    echo "[trimming] $SAMPLE_NAME"
    echo "java -jar $FOLDER/trimmomatic-0.39.jar PE ${ARR[$i]} ${ARR[$i+1]} -baseout ./trimmed/$SAMPLE_NAME.fastq.gz ILLUMINACLIP:$FOLDER/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:16"
    java -jar $FOLDER/trimmomatic-0.39.jar PE ${ARR[$i]} ${ARR[$i+1]} -baseout ./trimmed/$SAMPLE_NAME.fastq.gz ILLUMINACLIP:$FOLDER/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:16
}
cd trimmed
mkdir paired; mv *P.fastq.gz paired
mkdir unpaired; mv *U.fastq.gz unpaired
cd ..
date
