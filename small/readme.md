# Working with stranded NGS data
## Alignment to the reference with Bowtie2
### Paired-end reads
- BWA is unable to process stranded data (https://sourceforge.net/p/bio-bwa/mailman/message/31050505/) to go with only positive or negative strand data strand. The option is to use Bowtie for alignment.
- The following script will process multiple FASTQ.GZ files with Bowtie2. It should be two files with paired reads per sample (two FASTQ.GZ files) only given as input. No more, no less. Rename all these files to have a similar beginning of the filename ending with the “\_“ symbol. E.g for the 10A sample keep “10A” as the first symbols: 10A\_r1.FASTQ.GZ and 10A\_r2.FASTQ.GZ. 10A\_.FASTQ.GZ and 10A\_1.FASTQ.GZ are also fine. Any file like 10A.\_FASTQ.GZ or 10A.FASTQ.GZ will be wrong.
- Output will consist of two sorted BAM files. with positive-strand (XXXX-pos.bam) and negative-strand (XXXX-neg.bam) reads. Non-aligned reads will not be saved.
- Before starting, don't forget to modify reference filename (REFERENCE).

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

## Checking the strandedness direction with Kallisto.
It is possible to construct the library in a such way that reverse and forward reads are swapped. In order to check the dataset it is possible to run Kallisto and to see which direction maps better to the transcriptome.

Select only the first 0.5 million reads from each FASTQ file.
```
zcat WT3_1P.uniq.fastq.gz | head -n 500000 | gzip > WT3_1P.test.fastq.gz
zcat WT3_2P.uniq.fastq.gz | head -n 500000 | gzip > WT3_2P.test.fastq.gz
```

Create index for the human genome index ("transcriptome.idx") with Kallisto (described well in the RNA-seq pathway). Then run the follwoing commands.

```
~/bin/kallisto/kallisto quant -i transcriptome.idx  -o test.un WT3_1P.test.fastq.gz WT3_2P.test.fastq.gz 
~/bin/kallisto/kallisto quant -i transcriptome.idx  -o test.rf WT3_1P.test.fastq.gz WT3_2P.test.fastq.gz --rf-stranded
~/bin/kallisto/kallisto quant -i transcriptome.idx  -o test.fr WT3_1P.test.fastq.gz WT3_2P.test.fastq.gz --fr-stranded

paste test.fr/abundance.tsv test.rf/abundance.tsv test.un/abundance.tsv  | cut -f1,4,9,14 | awk 'BEGIN{sum1=0;sum2=0;sun3=0}{sum1+=$2;sum2+=$3;sum3+=$4}END{print sum1,sum2,sum3}' | less | awk '{print $2/$1,$3/$1,$3/$2}'| awk '{if($1<0.3 && $3>3)print "forward-stranded";else if($1>3 && $2>3)print "reverse-stranded";else print "unstranded"}'
```

# Basic Linux commands
## Creating MD5 sum for each individual file in the directory
```
find -type f -exec md5sum "{}" +
```
## Comressing all individual files in the directory into individual GZ archives
```
parallel gzip -v9 ::: *
```
## Checking number of aligned reads in BAM files with SeqKit
```
#!/bin/bash
TOTAL_FILES=`find -iname '*.bam' | wc -l`
ARR=($(ls *.bam))
for ((i=0; i<$TOTAL_FILES; i+=1)) {
    FILE_NAME=`echo ${ARR[$i]} | awk -F "." '{print $1}'`
    printf "\n"
    echo "[processing] $FILE_NAME"
    samtools index $FILE_NAME.bam
    seqkit bam -C $FILE_NAME.bam
}
```
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
The following script, in parallel, converts FASTQ into FASTA files and, therefore, deletes information on quality of the reads. Do quality trimming prior to running this script.
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

# Cleaning the Linux environment
- The following script will delete rubbish from the system. Use it with caution.
- Run it as a root.

```
#!/bin/sh

#Packages updating
sudo apt update
sudo apt-get upgrade

# Deleting partial packages
apt-get clean && apt-get autoclean
apt-get remove --purge -y software-properties-common

# Removing no longer required packages
apt-get autoremove -y

# Removing orphaned packages
deborphan | xargs sudo apt-get -y remove --purge

#deleting old kernels
sudo dpkg --list | egrep -i --color 'linux-image|linux-headers'
echo $(dpkg --list | grep linux-image | awk '{ print $2 }' | sort -V | sed -n '/'`uname -r`'/q;p') $(dpkg --list | grep linux-headers | awk '{ print $2 }' | sort -V | sed -n '/'"$(uname -r | sed "s/\([0-9.-]*\)-\([^0-9]\+\)/\1/")"'/q;p') | xargs sudo apt-get -y purge
sudo dpkg --list | egrep -i --color 'linux-image|linux-headers'

# Cleaning /tmp
find /tmp -type f -atime +2 -mtime +2  | xargs  /bin/rm -f &&
find /tmp -type d -mtime +2 -exec /bin/rm -rf '{}' \; &&
find /tmp -type l -ctime +2 | xargs /bin/rm -f &&
find -L /tmp -mtime +2 -print -exec rm -f {} \;

# Cleaning Chromium browser cache
rm -r ~/.cache/chromium
rm -r ~/.config/chromium/Default/File\ System

#Cleaning Chrome browser cache
rm -r /home/*/.cache/google-chrome/

#cleaning images thumbnails
rm -r /home/*/.cache/thumbnails

# Cleaning the Trash
rm -rf /home/*/.local/share/Trash/*/**
rm -rf /root/.local/share/Trash/*/**

# Cleaning old snap versions
snap list --all | while read snapname ver rev trk pub notes; do if [[ $notes = *disabled* ]]; then sudo snap remove "$snapname" --revision="$rev"; fi; done

# Clean all the log file. Delete all .gz and rotated file
find /var/log -type f -regex ".*\.gz$" | xargs rm -Rf
find /var/log -type f -regex ".*\.[0-9]$" | xargs rm -Rf
journalctl --vacuum-time=10d

echo "Cleaning is completed" 
```

# Microbiome analysis (STUB)
## Clean reads with Trimmomatic
This step is described for the iVAR pathway.
## Create 16S rRNA database for Kraken
```
DATABASE=silva
cd ~/soft/kraken2/
./kraken2-build --special silva --db $DATABASE
```
## Run Kraken
```
#!/bin/bash

DATABASE=silva
FOLDER=~/bin/kraken2

date
TOTAL_FILES=`find ./ -maxdepth 1 -iname '*.fastq.gz' | wc -l`
ARR=($(ls *.fastq.gz))
for ((i=0; i<$TOTAL_FILES; i+=2))
{
  SAMPLE_NAME=`echo ${ARR[$i]} | awk -F "_" '{print $1}'`
  printf "\n"
  echo "[processing] $SAMPLE_NAME"
  echo "${FOLDER}/kraken2 --db ${FOLDER}/${DATABASE} --threads 4 --output ${SAMPLE_NAME}_mapping.txt --report ${SAMPLE_NAME}_report.txt --classified-out ${SAMPLE_NAME}_pos#.fastq --unclassified-out ${SAMPLE_NAME}_neg#.fastq --gzip-compressed --paired ${SAMPLE_NAME}_1P.fastq.gz ${SAMPLE_NAME}_2P.fastq.gz"
  ${FOLDER}/kraken2 --db ${FOLDER}/${DATABASE} --threads 4 --output ${SAMPLE_NAME}_mapping.txt --report ${SAMPLE_NAME}_report.txt --classified-out ${SAMPLE_NAME}_pos#.fastq --unclassified-out ${SAMPLE_NAME}_neg#.fastq --gzip-compressed --paired ${SAMPLE_NAME}_1P.fastq.gz ${SAMPLE_NAME}_2P.fastq.gz
}
echo "[Sorting files to directories and compressing files]..."
mkdir 0.input_reads && mv *.gz $_
gzip -9 *.fastq
mkdir 1.sorted_reads && mv *.gz $_
mkdir 2.mapping && mv *_mapping.txt $_
mkdir 3.report && mv *_report.txt $_
date
echo "[FINISHED]"
```
## Transferring results to Krona
Use KrakenTools (https://ccb.jhu.edu/software/krakentools/index.shtml) to convert Kraken's output to Krona file format.

```
find ./ -maxdepth 1 -iname '*.txt' | parallel -j 4 --progress ~/bin/KrakenTools-master/kreport2krona.py -r {} -o {.}.krona_report.txt
find ./ -maxdepth 1 -iname '*.txt' | parallel -j 4 --progress ~/bin//Krona-master/KronaTools/scripts/ImportText.pl {} -o {.}.krona_report.html
```
Open HTML files with web-brower
## Pavian
Install RStudio.
Then install Pavian.
```
install.packages("remotes")
remotes::install_github("fbreitwieser/pavian")
pavian::runApp(port=5000)
```
Run Pavian in RStudio to visualise the results.
