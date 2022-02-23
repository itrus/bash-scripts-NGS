#!/bin/bash

# filenames for input and output files
CDNA=Homo_sapiens.GRCh38.cdna.all.fa.gz
GTF=Homo_sapiens.GRCh38.105.gtf.gz
INDEX=homo105.idx

date
KALLISTO=~/bin/kallisto
$KALLISTO/kallisto index -i $INDEX $CDNA
$KALLISTO/kallisto inspect -g $GTF $INDEX
date
