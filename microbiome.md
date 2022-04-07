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
Open HTML files with web browser
## Pavian
Install RStudio.
Then install Pavian.
```
install.packages("remotes")
remotes::install_github("fbreitwieser/pavian")
pavian::runApp(port=5000)
```
Run Pavian in RStudio to visualise the results.
