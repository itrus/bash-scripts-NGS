# Working with stranded NGS data
## Alignment to the reference with Bowtie2
### Paired-end reads
- BWA is unable to process stranded data and to go with only positive or negative strand data strand (https://sourceforge.net/p/bio-bwa/mailman/message/31050505/). The option is to use Bowtie for alignment.
- The following script will process multiple FASTQ.GZ files with Bowtie2. It should be only two files with paired reads per sample (two FASTQ.GZ files) given as input. No more, no less. Rename all these files to have a similar beginning of the filename ending with the “\_“ symbol. E.g. for the 10A sample keep “10A” as the first symbols: 10A\_r1.FASTQ.GZ and 10A\_r2.FASTQ.GZ. 10A\_.FASTQ.GZ and 10A\_1.FASTQ.GZ are also fine. Any file like 10A.\_FASTQ.GZ or 10A.FASTQ.GZ will be wrong.
- Output will consist of two sorted BAM files. with positive-strand (XXXX-pos.bam) and negative-strand (XXXX-neg.bam) reads. Non-aligned reads will not be saved.
- Before starting, don't forget to modify the reference filename (REFERENCE).

```
#!/bin/bash

REFERENCE=ZIKV-PR.fasta

TOTAL_FILES=`find -iname '*.fastq.gz' | wc -l`
ARR=($(ls *.fastq.gz))
bowtie2-build -f $REFERENCE index
for ((i=0; i<$TOTAL_FILES; i+=2)) {
    FILE_NAME=`echo ${ARR[$i]} | awk -F "." '{print $1}'`
    printf "\n"
    echo "[processing] $FILE_NAME"
    echo "[aligning positive reads]"
    bowtie2 --threads 8 --norc -x index -1 ${ARR[$i]} -2 ${ARR[$i+1]} | samtools view -b -F 4 -q 0 | samtools sort -o $FILE_NAME-pos.bam
    echo "[aligning negative reads]"
    bowtie2 --threads 8 --nofw -x index -1 ${ARR[$i]} -2 ${ARR[$i+1]} | samtools view -b -F 4 -q 0 | samtools sort -o $FILE_NAME-neg.bam
}
```
### Single/unpaired reads
```
#!/bin/bash

REFERENCE=ZIKV-PR.fasta

TOTAL_FILES=`find -iname '*.fastq.gz' | wc -l`
ARR=($(ls *.fastq.gz))
bowtie2-build -f $REFERENCE index
for ((i=0; i<$TOTAL_FILES; i+=1)) {
    FILE_NAME=`echo ${ARR[$i]} | awk -F "." '{print $1}'`
    printf "\n"
    echo "[processing] $FILE_NAME"
    echo "[aligning positive reads]"
    bowtie2 --threads 8 --norc -x index -U ${ARR[$i]} | samtools view -b -F 4 -q 0 | samtools sort -o $FILE_NAME-pos.bam
    echo "[aligning negative reads]"
    bowtie2 --threads 8 --nofw -x index -U ${ARR[$i]} | samtools view -b -F 4 -q 0 | samtools sort -o $FILE_NAME-neg.bam
}
```

## Checking the strandedness direction with Kallisto
It is possible to construct the library in such a way that reverse and forward reads are swapped. To check the dataset, it is possible to run Kallisto and to see which direction maps better to the transcriptome.

Select only the first 0.5 million reads from each FASTQ file.
```
zcat WT3_1P.uniq.fastq.gz | head -n 500000 | gzip > WT3_1P.test.fastq.gz
zcat WT3_2P.uniq.fastq.gz | head -n 500000 | gzip > WT3_2P.test.fastq.gz
```

Create an index for the human genome index ("transcriptome.idx") with Kallisto (described well in the RNA-seq pathway). Then run the following commands.

```
~/bin/kallisto/kallisto quant -i transcriptome.idx  -o test.un WT3_1P.test.fastq.gz WT3_2P.test.fastq.gz 
~/bin/kallisto/kallisto quant -i transcriptome.idx  -o test.rf WT3_1P.test.fastq.gz WT3_2P.test.fastq.gz --rf-stranded
~/bin/kallisto/kallisto quant -i transcriptome.idx  -o test.fr WT3_1P.test.fastq.gz WT3_2P.test.fastq.gz --fr-stranded

paste test.fr/abundance.tsv test.rf/abundance.tsv test.un/abundance.tsv  | cut -f1,4,9,14 | awk 'BEGIN{sum1=0;sum2=0;sun3=0}{sum1+=$2;sum2+=$3;sum3+=$4}END{print sum1,sum2,sum3}' | less | awk '{print $2/$1,$3/$1,$3/$2}'| awk '{if($1<0.3 && $3>3)print "forward-stranded";else if($1>3 && $2>3)print "reverse-stranded";else print "unstranded"}'
```

## Pipeline for checking multiple paired files at once
```
#!/bin/bash
date
KALLISTO=~/bin/kallisto
INDEX=homo108.idx
GTF=Homo_sapiens.GRCh38.108.gtf.gz
TOTAL_FILES=`find -iname '*.fastq.gz' | wc -l`
ARR=($(ls *.fastq.gz))
mkdir temp
for ((i=0; i<$TOTAL_FILES; i+=2))
{
   SAMPLE_NAME=`echo ${ARR[$i]} | awk -F "_" '{print $1}'`
   printf "\n"
   echo "[processing] $SAMPLE_NAME"
   zcat ${ARR[$i]} | head -n 400000 | gzip > temp/R1.test.fastq.gz
   zcat ${ARR[$i+1]} | head -n 400000 | gzip > temp/R2.test.fastq.gz
   echo "testing unstranded"
   $KALLISTO/kallisto quant -i $INDEX -g $GTF -t 4  -o test.un temp/R1.test.fastq.gz temp/R2.test.fastq.gz
   echo "testing RF-stranded"
   $KALLISTO/kallisto quant -i $INDEX -g $GTF -t 4  -o test.rf temp/R1.test.fastq.gz temp/R2.test.fastq.gz --rf-stranded
   echo "testing FR-stranded"
   $KALLISTO/kallisto quant -i $INDEX -g $GTF -t 4  -o test.fr temp/R1.test.fastq.gz temp/R2.test.fastq.gz --fr-stranded
   echo "RESULT:"
   paste test.fr/abundance.tsv test.rf/abundance.tsv test.un/abundance.tsv  | cut -f1,4,9,14 | awk 'BEGIN{sum1=0;sum2=0;sun3=0}{sum1+=$2;sum2+=$3;sum3+=$4}END{print sum1,sum2,sum3}' | less | awk '{print $2/$1,$3/$1,$3/$2}'| awk '{if($1<0.3 && $3>3)print "forward-stranded";else if($1>3 && $2>3)print "reverse-stranded";else print "unstranded"}'
}
date 
```
