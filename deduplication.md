## Deduplication of similar reads created by PCR
### Nubeam-dedup for paired-end reads
The following script will process multiple FASTQ.GZ files with Bowtie2. It should be two files with paired reads per sample (two FASTQ.GZ files) only given as input. No more, no less. Rename all these files to have a similar beginning of the filename ending with the “\_“ symbol. E.g for the 10A sample keep “10A” as the first symbols: 10A\_r1.FASTQ.GZ and 10A\_r2.FASTQ.GZ. 10A\_.FASTQ.GZ and 10A\_1.FASTQ.GZ are also fine. Any file like 10A.\_FASTQ.GZ or 10A.FASTQ.GZ will be wrong.
```
#!/bin/bash
date
FOLDER=/home/user/bin
TOTAL_FILES=`find -iname '*.fastq.gz' | wc -l`
ARR=($(ls *.fastq.gz))
for ((i=0; i<$TOTAL_FILES; i+=2))
{
    SAMPLE_NAME=`echo ${ARR[$i]} | awk -F "_" '{print $1}'`
    echo $SAMPLE_NAME
    echo "$FOLDER/nubeam-dedup --in1 ${ARR[$i]} --in2 ${ARR[$i+1]} --gz 9"
    $FOLDER/nubeam-dedup --in1 ${ARR[$i]} --in2 ${ARR[$i+1]} --gz 9
    echo "" 
}
date
```
### Using fastx_collapser for single-end reads.
The following script, in parallel, converts FASTQ into FASTA files and, therefore, deletes information on the quality of the reads. Do quality trimming before running this script.
```
#!/bin/bash
TOTAL_FILES=`find -iname '*.fastq' | wc -l`
ARR=($(ls *.fastq))
for ((i=0; i<$TOTAL_FILES; i+=1)) {
    FILE_NAME=`echo ${ARR[$i]} | awk -F "." '{print $1}'`
    printf "\n"
    echo "[processing] $FILE_NAME"
    sed -n '1~4s/^@/>/p;2~4p' $FILE_NAME.fastq > $FILE_NAME.fasta
    ~/bin/fastx/fastx_collapser -i $FILE_NAME.fasta -o $FILE_NAME-collapsed.fasta
}
```
