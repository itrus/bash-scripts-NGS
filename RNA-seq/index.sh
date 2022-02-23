#!/bin/bash
date
INDEX=homo105.idx
CDNA=Homo_sapiens.GRCh38.cdna.all.fa.gz
GTF=Homo_sapiens.GRCh38.105.gtf.gz
KALLISTO=~/bin/kallisto
$KALLISTO/kallisto index -i $INDEX $CDNA
$KALLISTO/kallisto inspect -g $GTF $INDEX
date
