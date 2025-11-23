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
#!/usr/bin/env bash
# deb12-system-cleanup.sh – safe monthly housekeeping for Debian 12
# version 2025-07-01
set -euo pipefail
IFS=$'\n\t'

log() { printf '[%s] %s\n' "$(date +'%F %T')" "$1"; }
free_before=$(df --output=avail -B1M / | awk 'NR==2{print $1}')

log 'APT update / full-upgrade …'
sudo apt-get update -y
sudo apt-get full-upgrade -y # pulls kernel meta-pkgs
sudo snap refresh

log 'Removing packages no longer required …'
sudo apt-get autoremove --purge -y                # prunes superseded kernels too
sudo apt-get autoclean -y && sudo apt-get clean

log 'Purging orphaned libraries …'
if command -v deborphan >/dev/null 2>&1; then
    deborphan --guess-all -z | xargs -0 -r sudo apt-get purge -y
fi

log 'deleting old kernels'
sudo dpkg --list | egrep -i --color 'linux-image|linux-headers'
echo $(dpkg --list | grep linux-image | awk '{ print $2 }' | sort -V | sed -n '/'`uname -r`'/q;p') $(dpkg --list | grep linux-headers | awk '{ print $2 }' | sort -V | sed -n '/'"$(uname -r | sed "s/\([0-9.-]*\)-\([^0-9]\+\)/\1/")"'/q;p') | xargs sudo apt-get -y purge
sudo dpkg --list | egrep -i --color 'linux-image|linux-headers'

log '/tmp: deleting items older than 3 d …'
sudo find /tmp -xdev -mindepth 1 -mtime +3 -print0 | sudo xargs -0 -r rm -rf --

log 'Docker: pruning unused layers / volumes …'
if command -v docker >/dev/null 2>&1; then
    sudo docker system prune -af --volumes        # typically 600 MB–3 GB
fi

log 'Flatpak: removing old revisions …'
if command -v flatpak >/dev/null 2>&1; then
    sudo flatpak uninstall --unused -y
fi

log 'Cleaning Chrome browser cache'
rm -r ~/.cache/google-chrome/

log 'Cleaning old snap versions'
snap list --all | while read snapname ver rev trk pub notes; do if [[ $notes = *disabled* ]]; then sudo snap remove "$snapname" --revision="$rev"; fi; done

log 'Vacuuming system journal → 14 d / 300 MB …'
sudo journalctl --vacuum-time=14d --vacuum-size=300M

log 'Deleting stale crash dumps & coredumps …'
sudo rm -rf /var/crash/*.crash /var/lib/systemd/coredump/* || true

log 'Clearing language-package caches …'
for d in ~/.cache/pip ~/.cache/pipenv ~/.cache/composer ~/.npm ~/.cache/_cacache; do
    [ -d "$d" ] && rm -rf "$d"
done

log 'Removing user thumbnails & trash …'
if [ -n "${SUDO_USER-}" ]; then
    sudo -u "$SUDO_USER" bash -c 'rm -rf ~/.cache/thumbnails/* ~/.local/share/Trash/{files,info}/*'
fi

free_after=$(df --output=avail -B1M / | awk 'NR==2{print $1}')
freed=$(( free_after - free_before ))        # MiB
log "Freed ${freed} MiB"

log 'Cleanup completed.'
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
