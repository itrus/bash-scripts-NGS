# Demultiplexing
The following script will use Cutadapt to remove 3' adaptor stored in 3adapter.fasta file from a batch of FASTQ files.

First create 3adapter.fasta file.
```
>3adapter
NAGATCGGAAGAGCACACGTCTG
```
Now you are ready to remove it with the following script.
```
#!/bin/bash
TOTAL_FILES=`find -iname '*.fastq' | wc -l`
ARR=($(ls *.fastq))
for ((i=0; i<$TOTAL_FILES; i+=1)) {
    FILE_NAME=`echo ${ARR[$i]} | awk -F "." '{print $1}'`
    echo "[processing] $FILE_NAME"
    printf "\n"
    printf "Reads, thousands: "
    cat $FILE_NAME.fastq | echo $((`wc -l`/4/1000))
    printf "\n"
    ~/bin/cutadapt -a file:3adapter.fasta -o $FILE_NAME.clipped3.fastq $FILE_NAME.fastq
    printf "\n"
    printf "Reads, thousands: "
    cat $FILE_NAME.clipped3.fastq | echo $((`wc -l`/4/1000))
    printf "\n"
}
```
Removing 5' adapter and demultiplexing.

First create 5adapter.fasta file with all adapter variants. Put first the longer ones.
```
>adapter1
NNNCACTAGCN
>adapter2
NNNGTGAGCN
>adapter3
NNNAGAGCN 
```
Now you are ready to remove it and to demultipolex FASTQ files with the following script.

Zero-length reads will not be preserved (**-m 1** parameter).
```
#!/bin/bash
TOTAL_FILES=`find -iname '*.fastq' | wc -l`
ARR=($(ls *.fastq))
for ((i=0; i<$TOTAL_FILES; i+=1)) {
    FILE_NAME=`echo ${ARR[$i]} | awk -F "." '{print $1}'`
    echo "[processing] $FILE_NAME"
    printf "\n"
    ~/bin/cutadapt -m 1 -g file:5adapter.fasta -o {name}_$FILE_NAME.fastq $FILE_NAME.fastq
}
```
