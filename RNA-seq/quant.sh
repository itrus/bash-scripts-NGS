#!/bin/bash

INDEX=homo105.idx
GTF=Homo_sapiens.GRCh38.105.gtf.gz
KALLISTO=~/bin/kallisto

date
TOTAL_FILES=`find -iname '*.fastq.gz' | wc -l`
ARR=($(ls *.fastq.gz))
for ((i=0; i<$TOTAL_FILES; i+=2))
{
   SAMPLE_NAME=`echo ${ARR[$i]} | awk -F "_" '{print $1}'`
   printf "\n"
   echo "[processing] $SAMPLE_NAME"
   echo "$KALLISTO/kallisto quant -i $INDEX -g $GTF -t 4 -o $SAMPLE_NAME ${ARR[$i]} ${ARR[$i+1]}"
   $KALLISTO/kallisto quant -i $INDEX -g $GTF -t 4 -o $SAMPLE_NAME ${ARR[$i]} ${ARR[$i+1]}
}
date