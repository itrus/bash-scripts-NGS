# Detection of novel DNA viruses
## 1 SUMMARY
This manual describes a bioinformatics pipeline used to detect viral genomes and to identify SNVs (single nucleotide variants) in Next-Generation Sequencing (NGS) data
In the first step, paired-end reads from illumina NGS data are preprocessed (Figure 1).
 
![de_novo](https://user-images.githubusercontent.com/9166776/155354021-2285d582-0eb7-4f0d-8214-06166ef99840.png)

Figure 1. Processing

After confirming the expected coverage and depth in the PhiX sequencing control, Kraken is used to remove host genome reads (negative selection) and to detect viral genomes present in leftover reads (positive selection). _De novo_ assembly is made with Trinity and SPAdes. Both of them have strong and week sides. Thus, it is better to compare the output of both of them to build the consensus genome. SNVs are detected with VarScan package.


## 2 PREREQUISITES
You will need a Linux machine for data preprocessing (preferrably Ubuntu-based one). Final analysis is done in the R environment (R-Studio) in Windows or Linux.
Run all commands in the console/terminal/shell (Show Applications button in the lowest left corner > Search for “Terminal” or “Konsole”). Configure your Linux shell. Go to the Terminal preferences > Unnamed profile > Scrolling > Remove the checkmark from “Limit scrollback to ...”.
For running each command, you need to enter the proper directory. For example, if you need to use the directory /home/ngs1/temporary. To enter it with terminal, execute the following command:
```
cd /home/ngs1/temporary
```
After running each command:
If several commands were copy-pasted, check that the last command was executed. If not, press Enter to launch the last command.
- Read all output to check if it was finished without error messages.
- Copy-paste the output of the terminal session to a proper file (log.txt) and save it in the working folder for the current procedure.
For each analysis step you need an empty working directory. Bring input files, then after successful processing delete input files, but leave all newly created files (results) and reference files (e.g. FASTA-files, etc.). This allows to go back and to repeat any step starting with any initial point. Besides this, it comes possible to investigate differences between files created under different parameters/conditions.
- Create subfolders required for storing results of each procedure.

## 3 PREPARE THE LINUX ENVIRONMENT
First install all the require dependencies:
```
sudo apt autoremove
sudo apt install parallel gawk gzip default-jre
```
During installation, machine will ask admin's password.

### Kraken installation
The complete manual for Kraken is available online (https://github.com/DerrickWood/kraken2/wiki/Manual). Example of using it for virus detection: https://virologyj.biomedcentral.com/articles/10.1186/s12985-018-1001-z.
1. Download and install Kraken2 (https://github.com/DerrickWood/kraken2/archive/master.zip)
```
unzip kraken2-master.zip
mkdir ~/bin/kraken2
cd kraken2-master
./install_kraken2.sh ~/bin/kraken2/
```
2. Install dustmaker to allow masking of senseless repeats for Kraken
```
sudo apt-get install ncbi-blast+
```
3. Create host genome database for Kraken. For porcine genome:
   - Visit the website: https://www.ncbi.nlm.nih.gov/datasets/genomes/?txid=9823. Select **pig - Sscrofa 11.1 > Download > Select Genomic sequence (FASTA) > Download.**
   - Unpack the downloaded ZIP file (0.7 Gb) to a new folder.
```
unzip ncbi_dataset.zip
```
   - Use unpacked folder with the porcine genome (ncbi_dataset/data/GCF_000003025.6) to create database containing pig genome (**Pig**). This will take some time and space (67 min, 33.3 Gb).
```
DATABASE=Pig
date
cd ~/bin/kraken2
./kraken2-build --download-taxonomy --db $DATABASE
find ~/Downloads/ncbi_dataset/data/GCF_000003025.6 -name '*.fna' -print0 | xargs -0 -I{} -n1 ./kraken2-build --add-to-library {} --db $DATABASE
./kraken2-build --build --db $DATABASE --threads 4
date
```
4. Create **Vir** database containing Viral genomes (33 min, 31.7 Gb)
```
DATABASE=Vir
date
cd ~/bin/kraken2
./kraken2-build --download-taxonomy --db $DATABASE
./kraken2-build --download-library viral --db $DATABASE
./kraken2-build --build --db $DATABASE --threads 4
date
```
### Krona installation
Download and install Krona to ~/bin directory (https://github.com/marbl/Krona/archive/master.zip). The last command will download taxonomy data for Krona. It is required to have compatibility with Kraken.
```
unzip -d ~/bin/ Krona-master.zip
mkdir ~/bin/Krona-master/KronaTools/taxonomy
sudo ~/bin/Krona-master/KronaTools/install.pl
~/bin/Krona-master/KronaTools/updateTaxonomy.sh
```
Install curl. Krona needs it to download taxonomy data.
```
sudo apt-get install curl
```
### Trimmomatic instalation
Install Trimmomatic (http://www.usadellab.org/cms/?page=trimmomatic). Download and unpack it to /home/user/bin/Trimmomatic-0.40/. UGENE contains a build-in version of Trimmomatic and can replace the Linux version.
```
unzip -d ~/bin/ ./Trimmomatic-0.40.zip
```

### SPAdes

Download and install SPAdes (http://cab.spbu.ru/files/release3.14.1/SPAdes-3.14.1-Linux.tar.gz). The following command will unpack it immediately to ~/bin directory.
```
tar -xzf SPAdes-3.14.1-Linux.tar.gz -C ~/bin/
```

### Trinity
It is important that Trinity and processed BAM files are placed at Linux partition. Both installation and files processing will fail if it any component is outside (e.g. on USB disk or external hard disk).

Install required packages.
```
sudo apt install jellyfish salmon bowtie2 cmake
```
Download and install Trinity (https://github.com/trinityrnaseq/trinityrnaseq/releases). The following command will unpack it directly to ~/bin directory.
```
tar -xzf trinityrnaseq-v2.11.0.FULL.tar.gz -C ~/bin/
```
Install it.
```
cd ~/bin/trinityrnaseq-v2.11.0/
sudo make
```

### R and RStudio
To install RStudio, download from the official website (https://rstudio.com/products/rstudio/download/#download) version for Ubuntu 18 (or your version). Then install the downloaded file (**Right mouse button > Open With Software Install > Install**)

Make sure you have the latest version of R. To manually reinstall R execute: 
```
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo apt update
sudo apt install r-base r-base-core r-recommended r-base-dev
```

To install RStudio, download from the official website (https://rstudio.com/products/rstudio/download/#download) version for Ubuntu 18 (or your version). Then install the downloaded file (**Right mouse button > Open With Software Install > Install**)

## 4 Detection of novel DNA viruses
## 4.1 Trimmomatic
Prepare original FASTQ.GZ files with Trimmomatic.

Copy-paste to a new working folder FASTQ.GZ files containing paired reads and /main/trimmomatic.sh script. Launch the script to proceed with Trimmomatic.
- This workflow is only valid for the paired-ended reads (PE). You cannot use it for single-ended (SE) reads.
- Default workflow parameters: **ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36**. Adapt it in the trimmomatic.sh script if necessary.
- Illimina HiSeq is using Phred33 quality scores and needs TruSeq3-PE-2.fa file from the Trimmomatic package for removing adaptors. Input files for adaptors may vary depending on the machine, etc. Check it each time.
- Phred >20 parameter used for quality trim corresponds to 99% accuracy in nucleotides.
- Adapt the directory pathway for the Trimmomatic package in the trimmomatic.sh script if necessary.
- After trimming, check the report of Trimmomatic carefully for errors and save it as a log file (log.txt). The output is 4 files sorted into two directories:
  - 2 files with paired reads (forward and reverse): \<SAMPLE\>\_1P.fastq.gz and \<SAMPLE\>\_2P.fastq.gz stored in ./trimmed/paired directory.
  - 2 files with unpaired reads (forward and reverse): \<SAMPLE\>\_1U.fastq.gz and \<SAMPLE\>\_2U.fastq.gz stored in ./trimmed/unpaired directory.
- Collect and keep generated output files from the newly created folder named ‘trimmed’. For the following steps, you will need only 2 files with paired reads.
	
### 4.2 Running Kraken for the first time (negative selection)
First we do negative selection for porcine genome. Kraken takes two files as input, thus it is enough to provide only the constant part of the filename (e.g. PCV3-A). As soon as it is negative selection you will need unsorted/negative fasta files to continue (PCV3-A_neg_1.fastq and PCV3-A_neg_2.fastq). If you bring files with other extension (e.g. FASTA or FA), then adjust it.
```
DATABASE=Pig
FILENAME=PCV3-A

~/bin/kraken2/kraken2 --db ~/soft/kraken2/$DATABASE --threads 4 --output ${FILENAME}_mapping.txt --report ${FILENAME}_report.txt --classified-out ${FILENAME}_pos#.fastq --unclassified-out ${FILENAME}_neg#.fastq --gzip-compressed --paired ${FILENAME}_1P.fastq.gz ${FILENAME}_2P.fastq.gz
```

Results for each processed file are
- useless mapping file (PCV3-A_mapping.txt)
- classification report file (PCV3-A_report.txt)
- positively selected fasta-files (PCV3-A_pos_1.fastq and PCV3-A_pos_2.fastq). Delete them. It is just pig genome.
- negatively selected fasta-files (PCV3-A_neg_1.fastq and PCV3-A_neg_2.fastq). Save them. We will need them to continue.

### 4.3 Running Kraken for the second time (positive selection)
Use negativele selected FASTQ files and viral database to make virus identification with Kraken.
```
DATABASE=Vir
FILENAME=PCV3-A_neg

~/bin/kraken2/kraken2 --db ~/soft/kraken2/kraken2/$DATABASE --threads 4 --output ${FILENAME}_mapping.txt --report ${FILENAME}_report.txt --classified-out ${FILENAME}_pos#.fastq --unclassified-out ${FILENAME}_neg#.fastq --paired ${FILENAME}_1.fastq ${FILENAME}_2.fastq
```

Results for each processed file are the same as above. But now classification file (PCV3-A_neg_report.txt) is useful.

### 4.4 Krona
Process classification report files provided by Kraken2 with Krona.

To process a single file:
```
~/bin/Krona-master/KronaTools/scripts/ImportTaxonomy.pl -m 3 -t 5 PCV3-A_neg_report.txt
```
To process all the TXT files in the current folder.
```
find ./ -maxdepth 1 -iname '*.txt' | parallel -j 4 --progress '~/soft/Krona-master/KronaTools/scripts/ImportTaxonomy.pl -m 3 -t 5  {} -o {.}.krona_report.html'
```

The report is in HTML format. Open it with web-browser.

## 5 Genome assembly

I've got better results with Trinity then with SPAdes.

### 5.1 SPAdes

Process files one-by-one. Change variable $SAMPLE for each sample you have. E.g. PCV3-B for PCV3-B_neg_pos_2.fastq and PCV3-B_neg_pos_2.fastq files.
```
SAMPLE=PCV3-B
~/bin/SPAdes-3.14.1-Linux/bin/spades.py --isolate -1 ${SAMPLE}_neg_pos_1.fastq -2 ${Sample}_neg_pos_2.fastq -o ${SAMPLE}
```
Resulting scaffolds could be found in **./SAMPLE-NAME/scaffolds.fasta file**. Use it find the correct sequence of your virus. Use reference from Genbank (accession numbers start with NC_) for assistance. E.g. NC_031753.1 for PCV3.

### 5.2 Trinity
It is important that Trinity and processed BAM files are placed at Linux partition. Both installation and files processing will fail if it any component is outside (e.g. on USB disk or external hard disk).

Run it.
```
SAMPLE=PCV3-B
~/bin/trinityrnaseq-v2.11.0/Trinity --seqType fq --max_memory 10G --left ${SAMPLE}_neg_pos_1.fastq --right ${SAMPLE}_neg_pos_2.fastq --CPU 6 --output trinity_${SAMPLE}
```
Resulting scaffolds could be found in **./trinity_SAMPLE-NAME/Trinity.fasta file**. Use it find the correct sequence of your virus. Use reference from Genbank (accession numbers start with NC_) for assistance. E.g. NC_031753.1 for PCV3.

## 6 Coverage analysis

Use FASTA file from the previous step to make index with bwa.
```
bwa index consensus_SK18.fas
```
Align FASTQ files to the reference. Better to use files after Trimmomatic. Coverage will be slightly (0.6-2.5%) higher
```
INDEX=consensus_SK18.fas

date
TOTAL_FILES=`find -name '*.fastq.gz' | wc -l`
ARR=($(ls *.fastq.gz))
for ((i=0; i<$TOTAL_FILES; i+=2))
{
printf "\n"
SAMPLE_NAME=`echo ${ARR[$i]} | awk -F "_" '{print $1}'`
echo "[mapping running for $[i/2+1]/$[TOTAL_FILES/2]] $SAMPLE_NAME"
echo "bwa mem -t 32 $INDEX ${ARR[$i]} ${ARR[$i+1]} | samtools view -b -F 4 -F 2048 | samtools sort -o $SAMPLE_NAME.bam"
bwa mem -t 32 $INDEX ${ARR[$i]} ${ARR[$i+1]} | samtools view -b -F 4 -F 2048 | samtools sort -o $SAMPLE_NAME.bam
}
date
```

Create the GENOME file. And then it will be used to create BED file with coverage. In the last step it will be converted to TXT file.
```
INDEX=consensus_SK18.fas

samtools faidx $INDEX
awk {'printf $1 "\t" $2'} $INDEX.fai > $INDEX.genome
find ./ -iname '*.bam' | parallel -j 4 --progress 'bedtools bamtobed -i {} > {.}.bed | echo {/}'
find ./ -iname '*.bed' | parallel -j 4 --progress 'bedtools genomecov -d -i {} -g consensus_SK18.fas.genome > {.}.txt'
```

## 7 Making phylogenetic tree
Download all full length genomes pf PCV3:
1. Go to: GenBank (https://www.ncbi.nlm.nih.gov/nuccore). 
2. Search for: (txid1868221\[Organism:noexp\] AND ( "1900"\[SLEN\] : "2200"\[SLEN\] ))
3. **Send to: >  Complete Record; Choose Destination: File; Format: FASTA; Sort by: Date Released > Create File**
4. Add your consensus sequence to the resulting FASTA file.
5. Align it. E.g. with MAFFT at the https://usegalaxy.org/ website.
   - Upload it (Download from URL or upload files from disk button)
   - Align it (MAFFT)
6. MEGA X
   - This is circular genome. If necessary transfer beginning to the end. To match sequences to the reference one.
   - Delete duplicated sequences
   - **Data > Phylogenetic Analysis > No**; **Data > Exit AlnExplorer**; **Distance > Compute pairwise distance**
   - **Export to CSV formatted file button > Export type: Column > Print/Save Matrix**
7. Import in Excel. Find pairs with distance of 0.0000% and delete one from each pair from the FASTA file. Rename strains. Open with WordPad. Search and delete similar parts of strain names (e.g. search and replace “PCV3” with “_”, “complete genome” with “_” and so on. Then at last search and replace “\_\_” with “\_”).
8. Delete (filter) low quality sequences (with huge insertions or deletions, not aligned well).
9. **Data > Phylogenetic Analysis**; **Data > Exit AlnExplorer**; **Phylogeny > Construct/Test NJ tree > Compute**

## 8 SNV detection

1. Download VarScan (https://github.com/dkoboldt/varscan).
2. Run it for each sample you have.
```
samtools mpileup -f consensus_SK18.fas PCV3-B.bam | java -jar VarScan.v2.4.4.source.jar pileup2snp
```

## 9 Errors to handle
### Kraken2 (downloading error: _rsync\_from_ncbi.pl: unexpected FTP path (new server?) for na_)
Error is described here: https://github.com/DerrickWood/kraken/issues/114
You need to add 3 lines to rsync_from_ncbi.pl

```
if ( $full_path =~/^na/){
next
}
```
The final file will look like this:

![kraken2_error](https://user-images.githubusercontent.com/9166776/155363735-de0e6f17-627a-4021-b9cf-6559077592b7.png)

### Kraken2 (problem with deleting processed files prior to database build-up)
You should not try building databases on remoted folders or USB sticks. Files should be on your local disc with Linux-native file system. Error message:
```
user@user-VirtualBox:~/Desktop/WinDownloads/kraken2$ ./kraken2-build --download-library bacteria --db $Database
Step 1/2: Performing rsync file transfer of requested files
Rsync file transfer complete.
Step 2/2: Assigning taxonomic IDs to sequences
Processed 22240 projects (49744 sequences, 90.87 Gbp)... done.
All files processed, cleaning up extra sequence files...rm: cannot remove 'all/GCF/000/003/215/GCF_000003215.1_ASM321v1/GCF_000003215.1_ASM321v1_genomic.fna.gz': Operation not permitted
rm: cannot remove 'all/GCF/000/003/645/GCF_000003645.1_ASM364v1/GCF_000003645.1_ASM364v1_genomic.fna.gz': Operation not permitted
…

…
rsync_from_ncbi.pl: can't clean up all/ directory: 256
```

### SPAdes: Resulting files from Kraken2 are automatically renamed to FASTQ.
This is important for SPAdes. It will crush if you bring the same file but with other file extension (e.g. *.fa or *.fasta). Error message:
```
user@user-VirtualBox:~/Desktop/WinDownloads$ ~/soft/SPAdes-3.14.1-Linux/bin/spades.py --isolate -1 PCV3-B_neg_pos_1.fa -2 PCV3-B_neg_pos_2.fa -o PCV3-B
Traceback (most recent call last):
  File "/home/user/soft/SPAdes-3.14.1-Linux/bin/spades.py", line 643, in <module>
    main(sys.argv)
  File "/home/user/soft/SPAdes-3.14.1-Linux/bin/spades.py", line 581, in main
    cfg, dataset_data, command_line = parse_args(args, log)
  File "/home/user/soft/SPAdes-3.14.1-Linux/bin/spades.py", line 238, in parse_args
    secondary_filling=False, restart_from=False)
  File "/home/user/soft/SPAdes-3.14.1-Linux/share/spades/spades_pipeline/options_parser.py", line 1072, in parse_args
    load_processed_dataset, restart_from, options)
  File "/home/user/soft/SPAdes-3.14.1-Linux/share/spades/spades_pipeline/options_parser.py", line 1012, in postprocessing
    support.check_dataset_reads(dataset_data, (args.only_assembler or args.rna), args.iontorrent, log)
  File "/home/user/soft/SPAdes-3.14.1-Linux/share/spades/spades_pipeline/support.py", line 689, in check_dataset_reads
    (key, id + 1, reads_library["type"]), log)
  File "/home/user/soft/SPAdes-3.14.1-Linux/share/spades/spades_pipeline/support.py", line 140, in check_file_not_empty
    if next(reads_iterator, None) is None:
  File "/home/user/soft/SPAdes-3.14.1-Linux/share/spades/spades_pipeline/common/SeqIO.py", line 108, in parse_fasta
    assert (rec_id[0] == '>')
AssertionError
```
