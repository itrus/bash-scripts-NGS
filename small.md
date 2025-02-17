# Basic Linux commands
## Processing all files in the directory
```
for i in *.fastq; do pyFastqDuplicateRemover.py -f $i -o collapsed_$i; done
```
## Creating MD5 sum for each file in the directory
```
find -type f -exec md5sum "{}" +
```
## Compressing all individual files in the directory into individual GZ archives
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
## Converting all BAM files to SAM files
```
parallel --plus 'samtools view -h {} -o {...}.sam' ::: *.bam
```
## Converting all BAM files to BAM files with either negative or positive strand reads
```
#!/bin/bash
mkdir pos-neg
find ./ -maxdepth 1 -iname '*.bam' | parallel -j 4 --progress 'samtools view -b -F 4 -F 2048 -F 16 {} -o ./pos-neg/{.}_pos.bam | echo {/}'
find ./ -maxdepth 1 -iname '*.bam' | parallel -j 4 --progress 'samtools view -b -F 4 -F 2048 -f 16 {} -o ./pos-neg/{.}_neg.bam | echo {/}'
```
# Cleaning the Linux environment
- The following script will delete rubbish from the system. Use it with caution.
- Run it as a root.

```
#!/bin/sh

# Packages updating
sudo apt update
sudo apt-get upgrade

# Deleting partial packages
apt-get clean && apt-get autoclean
apt-get remove --purge -y software-properties-common

# Removing no longer required packages
apt-get autoremove -y

# Removing orphaned packages
deborphan | xargs sudo apt-get -y remove --purge

# Deleting old kernels
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

# Cleaning Chrome browser cache
rm -r /home/*/.cache/google-chrome/

# Cleaning images thumbnails
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

## Writing console output to log file
```
./qc.sh 2>&1 | tee -a log.txt
```

## Merging paired FASTQ files

* this work is done without the use of _parallel_ package but still is using 12 threads
* input is 4 (2x2) paired FASTQ files per each sample
     - e.g. **SAMPLE1-readA_1P.fastq.gz, SAMPLE1-readA_2P.fastq.gz, SAMPLE1-readB_1P.fastq.gz, SAMPLE1-readB_1P.fastq.gz** for Sample #1 that was sequenced twice (readA and readB) and trimmed to have only paired reads
     - trimmining and merging here the paired reads is a mandatory preliminary step, because we need exactly the same number of reads in 1P and 2P files to keep them paired after merging readA and repeated readB
* output is 2 files per each sample
     - e.g. **SAMPLE1_readA_1P.fastq.gz.merged **and** SAMPLE1_readA_1P.fastq.gz.merged**
     - while the filename still contains "readA" it encompasses now "readA"+"readB" as shown with ".merged" in the filename. Consequnetly "readA" and ".merged" could be removed from the filenames of the newly created files.
* QC
     - feel free to run QC before and after merging to see that 3 paired reads (readA) + 2 paired reads (readB) = 5 paired reads for the merged file
     - feel free to run kallisto before and after merging to see that 2 pseudoaligned reads (readA) + 1 pseudoaligned read (readB) = 3 pseudoaligned reads for the merged file
    
```
#!/bin/bash
date

# Collect all .gz files into an array
ARR=($(ls *.gz))
TOTAL_FILES=${#ARR[@]}

# Set how many concurrent merges you want (e.g., 12 for a 12-core CPU)
max_jobs=12
running=0

# Loop through the files in chunks of 4
for ((i=0; i<$TOTAL_FILES; i+=4)); do
  {
    SAMPLE_NAME=$(echo "${ARR[$i]}" | awk -F "_" '{print $1}')
    echo -e "\n[merging] $SAMPLE_NAME"

    echo "zcat ${ARR[$i]} ${ARR[$i+2]} | gzip -v9 > ${ARR[$i]}.merged"
    zcat "${ARR[$i]}" "${ARR[$i+2]}" | gzip -v9 > "${ARR[$i]}.merged"

    echo "zcat ${ARR[$i+1]} ${ARR[$i+3]} | gzip -v9 > ${ARR[$i+1]}.merged"
    zcat "${ARR[$i+1]}" "${ARR[$i+3]}" | gzip -v9 > "${ARR[$i+1]}.merged"
  } &

  # Enforce max concurrency (12 jobs)
  ((running++))
  if [ $running -ge $max_jobs ]; then
    wait -n  # wait for exactly one job to finish
    ((running--))
  fi
done

wait  # wait for any remaining background jobs
echo -e "\nAll merging tasks complete."
```
