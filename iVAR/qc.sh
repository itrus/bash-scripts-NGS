#!/bin/bash
FOLDER=/home/user/bin/FastQC
SEARCH_PATTERN='*.fastq.gz'

TOTAL_FILES=`find ./ -maxdepth 1 -iname "*$SEARCH_PATTERN*" | wc -l`
ARR=($(ls $SEARCH_PATTERN))
echo "[running QC]"
ls $SEARCH_PATTERN | parallel -j 4 --progress $FOLDER/fastqc {}
printf "\n"
for ((i=0; i<$TOTAL_FILES; i+=1)) {
    FILE_NAME=`echo ${ARR[$i]} | awk -F "." '{print $1}'`    
    echo "[processing] $FILE_NAME"
    printf "Reads, thousands: "
    zcat ${ARR[$i]} | echo $((`wc -l`/4/1000)) 
} 
printf "\n"
printf "Total reads, thousands: "
zcat $SEARCH_PATTERN | echo $((`wc -l`/4/1000))
