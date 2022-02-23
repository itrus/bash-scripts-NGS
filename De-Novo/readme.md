# Detection of novel DNA viruses
## 1 SUMMARY
This manual describes a bioinformatics pipeline used to detect viral genomes and to identify SNVs (single nucleotide variants) in Next-Generation Sequencing (NGS) data
In the first step, paired-end reads from illumina NGS data are preprocessed (Figure 1).
 
Figure 1. Processing
After confirming the expected coverage and depth in the PhiX sequencing control, Kraken is used to remove host genome reads (negative selection) and to detect viral genomes present in leftover reads (positive selection). De novo assembly is made with Trinity and SPAdes. Both of them have strong and week sides. Thus it is better to compare the output of both of them to build the consensus genome. SNVs are detected with VarScan package.
2	PREREQUISITES
We use Linux for data preprocessing.
Run all commands in the console/terminal/shell (Show Applications button in the lowest left corner > Search for “Terminal” or “Konsole”). Configure your Linux shell. Go to the Terminal preferences > Unnamed profile > Scrolling > Remove the checkmark from “Limit scrollback to ...”.
For running each command, you need to enter the proper directory with the terminal session. For example, for trimming PhiX, you will need to use the directory /home/ngs1/ temporary. To enter it with terminal, execute the following command:
cd "/home/ngs1/temporary"
After running each command:
o	If several commands were copy-pasted, check that the last command was executed. If not, press Enter to launch the last command.
o	Read all output in Terminal to check if there are no errors.
o	Copy-paste the terminal session's output to a file (log.txt) and save it in the working folder for each step.
For each step, you need an empty working directory. Bring input files; after successful execution, delete input files, but preserve all newly created files (results) and reference files (e.g., FASTA-files, etc.). This allows us to go back and to repeat any step under the same or different conditions, if necessary.
Create subfolders for each analysis step:
Step	Subfolder name for each analysis step	Comments
Section # in this manual	Name		
na	Raw data	/0 raw data	
na	Working directory for temporary files	/temp	
4	Quality control check-up	/1 QC MiSeq	
5	Processing of the PhiX control	/2 QC PhiX	
5.1	Trimming PhiX	/2 QC PhiX/Trimmed	
5.2	Aligning PhiX to reference	/2 QC PhiX/BWA	
	Trimming	/3 Trimmed	
	Virus identification with Kraken	/4 Kraken	Subfolders contain reference for the host genome, results of negative and positive selections, and Krona visualisation report
	SPAdes	/5.1 SPAdes	
	Trinity	/5.2 Trinity	
	BAMs from BWA	/6 BWA	
	Coverage analysis	/7 Coverage	
	SNVs processing	/8 SNV	
	Phylogenetic tree analysis	/9 Tree	
3	PREPARE THE LINUX ENVIRONMENT
Kraken installation
The complete manual for Kraken is available online (https://github.com/DerrickWood/kraken2/wiki/Manual). Example of using it for virus detection: https://virologyj.biomedcentral.com/articles/10.1186/s12985-018-1001-z.
1.	Download and install Kraken2 (v 2.1.1; 0.2 Mb; https://github.com/DerrickWood/kraken2/archive/master.zip)
unzip kraken2-master.zip
mkdir ~/soft/kraken2
cd kraken2-master
./install_kraken2.sh ~/soft/kraken2/
2.	Install dustmaker to allow masking of senseless repeats for Kraken
sudo apt-get install ncbi-blast+
3.	Create host genome database for Kraken. For porcine genome:
a.	Visit the website: https://www.ncbi.nlm.nih.gov/datasets/genomes/?txid=9823
b.	Select pig - Sscrofa 11.1 > Download > Select Genomic sequence (FASTA) > Download
i.	678 Mb
c.	Unpack the downloaded file (678 Mb) to a new folder
unzip ncbi_dataset.zip
d.	Use unpacked folder with the porcine genome (ncbi_dataset/data/GCF_000003025.6) to create database containing pig genome (67 min, 33.3 Gb)
date
Database=Pig
cd ~/soft/kraken2
./kraken2-build --download-taxonomy --db $Database
find ~/Desktop/WinDownloads/ncbi_dataset/data/GCF_000003025.6 -name '*.fna' -print0 | xargs -0 -I{} -n1 ./kraken2-build --add-to-library {} --db $Database
./kraken2-build --build --db $Database --threads 4
date
4.	Create database containing Viral genomes (33 min, 31.7 Gb)
date
Database=Vir
cd ~/soft/kraken2
./kraken2-build --download-taxonomy --db $Database
./kraken2-build --download-library viral --db $Database
./kraken2-build --build --db $Database --threads 4
date
Krona installation
Download and install Krona to ~/soft directory (0.4 Mb; v.2.7.1; https://github.com/marbl/Krona/archive/master.zip). The last command will download taxonomy data for Krona. It is required to have compatibility with Kraken.
unzip -d ~/soft/ Krona-master.zip
mkdir ~/soft/Krona-master/KronaTools/taxonomy
sudo ~/soft/Krona-master/KronaTools/install.pl
~/soft/Krona-master/KronaTools/updateTaxonomy.sh

Install Trimmomatic (v 0.39; binary; http://www.usadellab.org/cms/?page=trimmomatic). Download and unpack it to ~/soft/Trimmomatic-0.39/. UGENE on Windows contains a build-in version of Trimmomatic and can replace a Linux version.
unzip ./Trimmomatic-0.39.zip
To install RStudio, download from the official website (https://rstudio.com/products/rstudio/download/#download) version for Ubuntu 18 (or your version) (https://download1.rstudio.org/desktop/bionic/amd64/rstudio-1.2.5033-amd64.deb). Then install the downloaded file (Right mouse button > Open With Software Install > Install)


---------------------------------------------

Detection of novel DNA viruses
5.	Prepare original FASTQ.GZ files with Trimmomatic
a.	use protocol from iVAR
6.	Run Kraken twice
a.	First, do negative selection for porcine genome. Kraken takes two files as input, thus it is enough to provide only the constant part of the filename (e.g. PCV3-A). As soon as it is negative selection you will need unsorted/negative fasta files to continue (PCV3-A_neg_1.fastq and PCV3-A_neg_2.fastq). File extensions for input and output files are highlighted with green. If you bring files with other extension (e.g. FASTA or FA), then adjust it.
Database=Pig
Filename=PCV3-A 
~/soft/kraken2/kraken2 --db ~/soft/kraken2/$Database --threads 4 --output ${Filename}_mapping.txt --report ${Filename}_report.txt --classified-out ${Filename}_pos#.fastq --unclassified-out ${Filename}_neg#.fastq --gzip-compressed --paired ${Filename}_1P.fastq.gz ${Filename}_2P.fastq.gz 
b.	Results for each processed file are
i.	useless mapping file (PCV3-A_mapping.txt)
ii.	classification report file (PCV3-A_report.txt)
iii.	positively sorted fasta-files (PCV3-A_pos_1.fastq and PCV3-A_pos_2.fastq). Delete them. It is just pig genome.
iv.	negatively sorted fasta-files (PCV3-A_neg_1.fastq and PCV3-A_neg_2.fastq). Save them. We will need them to continue.
c.	Use negative fasta-files and viral database to make virus identification with Kraken
Database=Vir
Filename=PCV3-A_neg 
~/soft/kraken2/kraken2 --db ~/soft/kraken2/kraken2/$Database --threads 4 --output ${Filename}_mapping.txt --report ${Filename}_report.txt --classified-out ${Filename}_pos#.fastq --unclassified-out ${Filename}_neg#.fastq --paired ${Filename}_1.fastq ${Filename}_2.fastq 
d.	Results for each processed file are the same as above. But now classification file (PCV3-A_neg_report.txt) is useful.
7.	Install Krona
a.	Install curl. Krona needs it to download taxonomy data.
sudo apt-get install curl
8.	Process classification report files provided by Kraken2 with Krona
a.	Single file
~/soft/Krona-master/KronaTools/scripts/ImportTaxonomy.pl -m 3 -t 5 PCV3-A_neg_report.txt
b.	All txt-files in the current folder
find ./ -maxdepth 1 -iname '*.txt' | parallel -j 4 --progress '~/soft/Krona-master/KronaTools/scripts/ImportTaxonomy.pl -m 3 -t 5  {} -o {.}.krona_report.html'
c.	The report is in HTML format. Open it with web-browser.

iVAR for genome assembly
Protocol:
1.	Download Fasta file for PCV3 reference from GenBank (NC_031753.1)
2.	Trimmomatic for FASTQ.GZ > BWA to the reference genome and sorting BAM > Coverage > Create the GFF file for the reference genome
a.	Highlighted with red steps are described for ZIKV.
3.	Run modified iVAR script (primers masking is excluded for this run)
#!/bin/bash
clear
date
REFERENCE=pcv3_ref.fas
GFF_FILE=pcv3_ref.gff3
TOTAL_FILES=`find -iname '*.bam' | wc -l`
ARR=($(ls *.bam))
bwa index $REFERENCE
for ((i=0; i<$TOTAL_FILES; i+=1))
    {
    SAMPLE_NAME=`echo ${ARR[$i]} | awk -F "." '{print $1}'`
    printf "\n"
    echo "[processing] $SAMPLE_NAME"
    echo "[indexing] samtools index $SAMPLE_NAME.bam"
    samtools index $SAMPLE_NAME.bam
    echo "[making consensus] samtools mpileup -A -d 1000000 -Q 0 $SAMPLE_NAME.bam | ivar consensus -p $SAMPLE_NAME.consensus"
    samtools mpileup -A -d 1000000 -Q 0 $SAMPLE_NAME.bam | ivar consensus -p $SAMPLE_NAME.consensus
    echo "[generating final TSV file] samtools mpileup -A -d 1000000 --reference $REFERENCE -Q 0 $SAMPLE_NAME.bam | ivar variants -p $SAMPLE_NAME.final -t 0.03 -r $REFERENCE -g $SGFF_FILE"
    samtools mpileup -A -d 1000000 --reference $REFERENCE -Q 0 $SAMPLE_NAME.bam | ivar variants -p $SAMPLE_NAME.final -t 0.03 -r $REFERENCE -g $GFF_FILE
    }
echo "[Sorting files to directories and renaming]"
mkdir 0.input && mv $PRIMERS $REFERENCE $PRIMER_PAIRS $GFF_FILE $_
mkdir 1.index_files && mv *.pac $_ && mv *.sa $_ && mv *.ann $_ && mv *.amb $_ && mv *.fai $_ && mv *.bai $_ && mv *.bwt $_
mkdir 2.Consensus && mv *.consensus.fa $_
mkdir 3.Quality && mv *.qual.txt $_
mkdir 4.Final-TSV && mv *.final.tsv $_
mkdir 0.BAM-original && mv *.bam $_
cd ./4.Final-TSV
for FILE in *; do mv "$FILE" "${FILE/final./}"; done
cd ..
date
echo "[FINISHED]"
4.	Now you’ve generated consensus files for novel PCV3 strain: 2. Consensus/PCV3-A.consensus. As soon as this file is the same for all samples (no difference) you can use just one to continue.
5.	BWA of FASTQ.GZ files after Trimmomatic to the consensus genome and sorting BAM > Coverage
6.	Run iVAR for the second time. Code is the same as above. Just replace reference genome (pcv3_ref.fa) with consensus sequence (PCV3-B.consensus.fa).
7.	Combine TSV files with variants to Excel file and analyse SNVs.

SPAdes
•	Manual: https://github.com/ablab/spades/blob/spades_3.14.1/README.md
Protocol:
1.	Download and install SPAdes (27 Mb; v 3.14.1; http://cab.spbu.ru/files/release3.14.1/SPAdes-3.14.1-Linux.tar.gz)
a.	This command will unpack it immediately to ~/soft directory.
tar -xzf SPAdes-3.14.1-Linux.tar.gz -C ~/soft/
2.	Process files one-by-one. Change variable $Sample for each sample you have. E.g. PCV3-B for PCV3-B_neg_pos_2.fastq and PCV3-B_neg_pos_2.fastq files.
Sample=PCV3-B
~/soft/SPAdes-3.14.1-Linux/bin/spades.py --isolate -1 ${Sample}_neg_pos_1.fastq -2 ${Sample}_neg_pos_2.fastq -o ${Sample}
3.	Resulting scaffolds could be found in ./SAMPLE-NAME/scaffolds.fasta file. Use it find the correct sequence of your virus. Use reference from Genbank (accession numbers start with NC_) for assistance. E.g. NC_031753.1 for PCV3.

Trinity
•	Manual: https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-Trinity
Protocol:
1.	It is important that Trinity and processed BAM files are placed at Linux partition. Both installation and files processing will fail if it any component is outside.
2.	Install required packages.
sudo apt install jellyfish salmon bowtie2 cmake
3.	Download and install Trinity (61 Mb; v 2.11.0; https://github.com/trinityrnaseq/trinityrnaseq/releases)
a.	This command will unpack it immediately to ~/soft directory.
tar -xzf trinityrnaseq-v2.11.0.FULL.tar.gz -C ~/soft/
a.	Install it.
cd ~/soft/trinityrnaseq-v2.11.0/
sudo make
4.	Run it
Sample=PCV3-B
~/soft/trinityrnaseq-v2.11.0/Trinity --seqType fq --max_memory 10G --left ${Sample}_neg_pos_1.fastq --right ${Sample}_neg_pos_2.fastq --CPU 6 --output trinity_${Sample}
5.	Resulting scaffolds could be found in ./trinity_SAMPLE-NAME/Trinity.fasta file. Use it find the correct sequence of your virus. Use reference from Genbank (accession numbers start with NC_) for assistance. E.g. NC_031753.1 for PCV3.
Coverage (BWA + BEDtools + R)
bwa index consensus_SK18.fas

date
INDEX=consensus_SK18.fas
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
6.	Better to use files after Trimmomatic. Coverage will be higher

Sample	Read depth for FASTA after Trimmomatic, x	Read depth for FASTA after Kraken2, x
PCV3-A	76.93	74.93
PCV3-B	640.05	636.64
PCV3-C	40.29	39.92
PF-NTC	0.00	0.00

INDEX=consensus_SK18.fas
samtools faidx $INDEX
awk {'printf $1 "\t" $2'} $INDEX.fai > $INDEX.genome

find ./ -iname '*.bam' | parallel -j 4 --progress 'bedtools bamtobed -i {} > {.}.bed | echo {/}'
find ./ -iname '*.bed' | parallel -j 4 --progress 'bedtools genomecov -d -i {} -g consensus_SK18.fas.genome > {.}.txt'

Phylogenetic tree
1.	Download all full length genomes pf PCV3
a.	Go to: GenBank (https://www.ncbi.nlm.nih.gov/nuccore)
b.	Search for: (txid1868221[Organism:noexp] AND ( "1900"[SLEN] : "2200"[SLEN] ))
c.	Send to: >  Complete Record; Choose Destination: File; Format: FASTA; Sort by: Date Released > Create File
d.	Resulting FASTA file (1.3 Mb, 647 genomes on 201208)
2.	Add your consensus sequence (total is 648 genomes).
3.	Align it.
a.	https://usegalaxy.org/
b.	Upload it (Download from URL or upload files from disk button)
c.	Align it (MAFFT)
4.	MEGA X
a.	This is circular genome. If necessary transfer beginning to the end. To match sequences to the reference one.
b.	Delete duplicated sequences
i.	Data > Phylogenetic Analysis > No
ii.	Data > Exit AlnExplorer
iii.	Distance > Compute pairwise distance
iv.	Export to CSV formatted file button > Export type: Column > Print/Save Matrix
v.	Import in Excel. Find pairs with distance of 0.0000% and delete one from each pair from the FASTA file.
c.	Rename strains. Open with WordPad. Search and delete similar parts of strain names (e.g. search and replace “PCV3” with “_”, “complete genome” with “_” and so on. Then at last search and replace “__” with “_”).
d.	Delete (filter) low quality sequences (with huge insertions or deletions, not aligned well). Final number of sequences: 642.
e.	Data > Phylogenetic Analysis
f.	Data > Exit AlnExplorer
g.	Phylogeny > Construct/Test NJ tree > Compute
5.	
SNV detection
find ./ -maxdepth 1 -iname '*.bam' | parallel -j 4 --progress 'bcftools mpileup -Ou -f consensus_SK18.fas {} | bcftools call --ploidy 1 -mv -Ov -o {.}.vcf'
1.	Download VarScan (v. 2.4.4; https://github.com/dkoboldt/varscan; 0.2 Mb).
2.	Run it for each sample you have.
samtools mpileup -f consensus_SK18.fas PCV3-B.bam | java -jar VarScan.v2.4.4.source.jar pileup2snp
Errors to handle
1.	Kraken2: Downloading error: rsync_from_ncbi.pl: unexpected FTP path (new server?) for na 
a.	Error is described here: https://github.com/DerrickWood/kraken/issues/114
b.	You need to add 3 lines to rsync_from_ncbi.pl
if ( $full_path =~/^na/){
next
}
c.	 The final file will look like this:
 
2.	Kraken2: Problem with deleting processed files prior to database build-up. You should not trying building databases on remoted folders or USB sticks. Files should be on your local disc with Linux-native file system.
a.	Error message
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

9.	SPAdes: Resulting files from Kraken2 are automatically renamed to FASTQ. This is important for SPAdes. It will crush if you bring the same file but with other file extension (e.g. *.fa or *.fasta). Error message:
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


---------------------------------------------


First install all require dependencies:
sudo apt autoremove
sudo apt install autotools-dev gcc zlib1g-dev libbz2-dev \
liblzma-dev make libncurses5-dev autoconf g++ bwa r-base \
bedtools parallel default-jre mc krusader krename kate gawk

To install HTSlib libraries, download HTSlib 1.10 (1 Mb; http://www.htslib.org/download/), then finish the installation procedure.
tar -xvf ./htslib-1.10.tar.bz2
cd ./htslib-1.10/
./configure
make
sudo make install
If HTSlib is installed in a non-standard location (e.g.,/usr/local/lib/) and iVAR cannot find it, please add the following line to your BASH configuration file (~/.bash_profile or ~/.bashrc) so that iVar can find HTSlib dynamic libraries during runtime. No restart is needed.
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/
To install SAMlib libraries, download SAMtools 1.10 (4 Mb; http://www.htslib.org/download/) then finish the installation procedure.
tar -xvf ./samtools-1.10.tar.bz2
cd ./samtools-1.10/
./configure
make
sudo make install
Install iVAR (https://github.com/andersen-lab/ivar). Download iVAR 20200605 (1 Mb; https://github.com/andersen-lab/ivar) then finish the installation procedure.
unzip ./ivar-master.zip
cd ./ivar-master/
./autogen.sh
./configure
make
sudo make install
4	QUALITY CONTROL
Use Illumina Seq Analysis Viewer 2.4.7.0 on Windows to see the general report of the sequencing (Browse > Use MiSeq output whole folder as input for Illumina Seq Analysis Viewer). The full manual for the software package is available at illumina’s website (https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/sav/sequencing-analysis-viewer-user-guide-15020619-f.pdf).
Create a report in Word containing charts and data. To enlarge the chart click the chevron in the upper right corner of each chart. To copy-paste charts, click it with the right mouse button and select Copy To Clipboard.
1.	Analysis tab
a.	Data by Cycle
i.	Error Rate (percentage of wrong bases in the PhiX control reads).
ii.	Called Int  (corrected intensity of base calls)
b.	QScore Distribution (Q20 corresponds to 1% error rate, Q30 corresponds to 0.1% error rate)
c.	Data by Lane
i.	% Aligned (Frequency of reads representing PhiX)
d.	Qscore Heatmap (Q-scores by cycle)
2.	Summary tab. Press Copy To Clipboard to copy-paste data.
a.	Total – Yield Total (G) (Total number of Gbp saved)
b.	Total – Aligned (%) (Frequency of reads representing PhiX)
c.	Error Rate (%) (percentage of wrong bases in the PhiX control reads)
3.	Indexing tab.
a.	Total reads (total number of reads saved)
b.	% reads Identified (PF) (percentage of reads with barcodes)
c.	Min and Max (percentage of the least and most frequent barcodes)
d.	Click on a table with individual indexes with left mouse button, then press Control+A > Control+C to copy-paste the table content.
e.	Barchart contains a graphical representation of the table. Copy-paste chart by clicking it with the right mouse button and selecting Copy To Clipboard.
Use the corresponding folder to save the results (screenshots).
5	PROCESSING OF THE PHIX CONTROL
5.1	Trim adaptors and do quality trim with Trimmomatic
From the illumina machine output folder, copy-paste to a new working folder 2 PhiX FASTQ.GZ files containing non-barcoded paired reads that should include PhiX.
•	./0 raw/191125_M04229_0101_000000000-CRM9R/Data/Intensities/BaseCalls/Undetermined_S0_L001_R1_001.fastq.gz
•	./0 raw/191125_M04229_0101_000000000-CRM9R/Data/Intensities/BaseCalls/Undetermined_S0_L001_R2_001.fastq.gz.
Run Trimmomatic.
o	This workflow is only valid for the paired-ended reads (PE). You cannot use it for single-ended (SE) reads.
o	Workflow parameters: ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
o	Illimina HiSeq is using Phred33 quality scores and needs TruSeq3-PE-2.fa file from the Trimmomatic package for removing adaptors. Input file for adaptors may vary depending on machine, etc. Check it each time.
o	Phred >20 parameter used for quality trim corresponds to 99% accuracy in nucleotides.
o	Adapt the directory pathway for the Trimmomatic package (highlighted) if necessary.
clear
date
mkdir trimmed
FOLDER=/home/ngs1/soft/Trimmomatic-0.39
TOTAL_FILES=`find -iname '*.fastq.gz' | wc -l`
ARR=($(ls *.fastq.gz))
for ((i=0; i<$TOTAL_FILES; i+=2))
{
SAMPLE_NAME=`echo ${ARR[$i]} | awk -F "_" '{print $1}'`
printf "\n"
echo "[trimming] $SAMPLE_NAME"
echo "java -jar $FOLDER/trimmomatic-0.39.jar PE ${ARR[$i]} ${ARR[$i+1]} -baseout ./trimmed/$SAMPLE_NAME.fastq.gz ILLUMINACLIP:$FOLDER/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36"
java -jar $FOLDER/trimmomatic-0.39.jar PE ${ARR[$i]} ${ARR[$i+1]} -baseout ./trimmed/$SAMPLE_NAME.fastq.gz ILLUMINACLIP:$FOLDER/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
}
cd trimmed
mkdir paired; mv *P.fastq.gz paired
mkdir unpaired; mv *U.fastq.gz unpaired
cd ..
date
After trimming, check the report of Trimmomatic carefully for errors and save it as a log file (log.txt). The output is 4 files sorted in two directories:
o	2 files with paired reads (forward and reverse): <SAMPLE>_1P.fastq.gz and <SAMPLE>_2P.fastq.gz stored in ./trimmed/paired directory.
o	2 files with unpaired reads (forward and reverse): <SAMPLE>_1U.fastq.gz and <SAMPLE>_2U.fastq.gz stored in ./trimmed/unpaired directory.
Collect and keep generated output files from the newly created folder named ‘trimmed’. For the following steps, you will need only 2 files with paired reads.
5.2	Align reads to the reference
Make a copy of output files from Trimmomatic corresponding to paired reads (2 files per sample). For the PhiX control sample name is Undetermined.
Provide a file representing the reference for the alignment.
•	This is an important step! Get the correct reference (FASTA file) for PhiX from the illumina website (https://support.illumina.com/sequencing/sequencing_software/igenome.html). Genbank has the PhiX sequence (NC_001422.1) that is 5 nt different from the illumina’s one.
Place the reference FASTA file in the same folder as FASTQ.GZ files from Trimmomatic. Then make index files for the reference genome (highlighted).
bwa index illumina-phix.fa
To process paired reads, use the following commands:
•	When using another reference genome, change index name (highlighted in yellow).
date
INDEX=illumina-phix.fa
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
The resulting output is a sorted and merged BAM file named as <SAMPLE>. BAM. Unmapped reads are skipped (“F” parameter). After successfully finalizing the run, delete the copy-pasted input FASTQ.GZ files to save space, but keep the generated output BAM file and logging information from the shell.
5.3	Coverage
To visually inspect the coverage, open the resulting BAM-file with UGENE.
PhiX genome should be covered fully (100% genome coverage) with a high coverage depth (>10’000x).
6	SNV ANALYSIS IN ZIKV GENOME
6.1	Trim adaptors and do quality trim with Trimmomatic
Do it exactly as for PhiX, but use FASTQ.GZ files corresponding to the experimental samples from the illumina machine output folder. Similarly to PhiX processing, after successful trimming, you will need only 2 files with paired reads for the following steps.
6.2	Align reads to the reference
This step is done exactly as for PhiX analysis. Make a copy of output files from Trimmomatic corresponding to paired reads (2 files per sample). The only difference is the reference genome provided. Similarly to PhiX processing, don’t mix paired and unpaired files after Trimmomatic. Proceed only with paired ones. Don’t forget to save the log file from Terminal.
•	This is an important step! Don’t forget to provide the CORRECT reference genome. Get the correct reference (FASTA file) (e.g., for the PR strain of ZIKV from the Genbank (KU501215.1)).
The resulting output is a sorted and merged BAM file named as <SAMPLE>. BAM.
6.3	Coverage
To visually inspect the coverage of any sample, open the resulting BAM-file with UGENE.
For batch processing of all samples at once and detailed analysis:
1.	 Create a GENOME file from the reference FASTA file (e.g., zika-pr.fasta). 
INDEX=zika-pr.fasta
samtools faidx $INDEX
awk {'printf $1 "\t" $2'} $INDEX.fai > $INDEX.genome
2.	Copy the final (merged and sorted) BAM-files to a new directory. Place it in the same folder where the correct GENOME file (highlighted in yellow). The following commands will produce TXT files with the coverage. Besides this BED files are created. You may delete them as soon as they will not be required further.
find ./ -iname '*.bam' | parallel -j 4 --progress 'bedtools bamtobed -i {} > {.}.bed | echo {/}'
find ./ -iname '*.bed' | parallel -j 4 --progress 'bedtools genomecov -d -i {} -g zika-pr.fasta.genome > {.}.txt'
3.	Analyse resulting files with R script.
a.	Create the R script (coverage.Rmd) containing following code.
---
title: "Coverage parsing"
author: "Ivan Trus"
date: "2020-10-19"
output:
  html_document: default
  pdf_document: default
---
Execution start time: *`r Sys.time()`*
```{r setup, include=FALSE}
knitr::opts_chunk$set(
  echo = FALSE,
  message = TRUE,
  warning = TRUE
)
library(vroom)
library(reshape)
```

# Initialization of the script and data loading
```{r init, message=FALSE}
Filelist <- list.files(pattern = ".*.txt")
CountGenomes <- length(Filelist)
Threshold <- 400

# Defining functions
ListGaps <- function(Gaps) {
  Start <- append(Gaps[1], Gaps[which(diff(Gaps) > 1) + 1])
  Start
  End <- append(Gaps[which(diff(Gaps) > 1)], Gaps[length(Gaps)])
  End
  GapsSummary <- data.frame(Start, End)
  colnames(GapsSummary) <- c("Start", "End")
  return(GapsSummary)
}

Coverage <- vroom(Filelist, id = Filelist, col_names = "")
Coverage <- Coverage[c(1, 3, 4)]
colnames(Coverage)[1] <- "X1"
Coverage <- t(cast(Coverage, X1 ~ X2, value.var = "X3"))
colnames(Coverage) <- gsub(".txt", "", Filelist)
GenomeLength <- length(Coverage[, 1])

```

* Current working directory: **`r getwd()`**
* List of files to process: **`r Filelist`**
* Total number of files: **`r CountGenomes`**
* Threshold for coverage counting: **`r Threshold`** x
* Genome length: **`r GenomeLength`** nt

# Diagnostic data, chart and list of gaps for each sample
```{r main}
Results <- matrix(
  nrow = CountGenomes, ncol = 3,
  dimnames = list(1:CountGenomes, c(
    "Sample", "Read depth, x", "Coverage, %"))
)
Results <- as.data.frame(Results)
Results[, 1] <- colnames(Coverage)
AllGaps <- NA

summary(Coverage)

for (i in 1:CountGenomes) {
  Results[i, 2] <- round(mean(Coverage[, i]), digits = 4)
  Results[i, 3] <- round(100 - sum(Coverage[, i] < Threshold) / GenomeLength * 100, digits = 4)

  # We need to replace 0 with 1 to make it possible to be plotted fully on log scale
  GapsForLogChart <- Coverage[, i]
  GapsForLogChart[GapsForLogChart == 0] <- 1

  if(max(GapsForLogChart) < Threshold) YaxisLimit <- Threshold else YaxisLimit <- max(GapsForLogChart)
  plot(log10(GapsForLogChart),
    type = "l", xlab = "Nucleotide position in the ZIKV genome",
    ylab = "Read depth (log10)", ylim = c(0, log10(YaxisLimit) * 1.05))
  text(x = 0, y = 1.05 * log10(YaxisLimit), colnames(Coverage)[i], pos = 4, font = 2, cex = 1.2)
  abline(h = log10(Threshold), col = "red", lwd = 3)
  GapsPositions <- as.numeric(names(subset(Coverage[, i], Coverage[, i] < Threshold)))
  AllGaps <- sort(c(AllGaps, setdiff(GapsPositions, AllGaps)))
  print(paste("The list of gaps in", colnames(Coverage)[i], "with the depth less than", Threshold, "x"))
  print(ListGaps(GapsPositions))
}
```

# Summary on read depth and coverage
```{r final}
Results
print(paste("Combined gaps with coverage depth less than", Threshold, "x."))
print(ListGaps(AllGaps))
print(paste("Intersection of coverage in all provided files:", round(100 - length(AllGaps) / GenomeLength * 100, digits = 4), "%"))
```

Execution end time: *`r Sys.time()`*
b.	X and Y graph axis titles are highlighted and can be adapted as needed.
c.	Place the created R script (coverage.Rmd) and TXT files with coverages to a separate (new) directory. Make sure that only the relevant TXT files were placed in this folder.
d.	Open and run the R script with RStudio.
i.	Use Code > Run Region > Run All (Control+Alt+R) if you want to see results in RStudio.
ii.	Use File > Knit Document if you need a report (HTML file) generated. This file can be opened with any web browser. To create a report in PDF file format, you need first to install the LaTeX distribution. 
iii.	If there is an Error: “Line 22 Error in loadnamesspace(I, c(lib.loc, .libPaths()), versioncheck = vI[[i]]) : namespace ‘rlang’ 0.4.6 is already loaded, but >= 0.4.7 is required calls: <Anonymous> … getNamespace -> namespaceImport -> loadNamespace”, update ‘rlang’ package in R studio.
6.4	Merge files representing one sample but pool #1 and pool #2
Do this step when DNA representing the whole ZIKV or PCV2 genome (amplified with two different primer pools; see materials and methods for PrimalSeq NGS) were not mixed during NGS library preparation. If pools were mixed, each biological sample is represented by one BAM file only, and merging is not required. Pools are not the same as technical replicates.
1.	Rename all files manually. E.g. for Pool #1 and Pool #2: files 4A.bam and 5A.bam should be renamed to 4A.bam and 4A_5A.bam.
2.	The following commands process pair of files at once. Thus, temporarily remove singlets (negative control samples represented by pool #1 or pool #2 only) to a separate folder.
mkdir 1-merged
mkdir 2-sorted
TOTAL_FILES=`find -iname '*.bam' | wc -l`
ARR=($(ls *.bam))
for ((i=0; i<$TOTAL_FILES; i+=2))
{
SAMPLE_NAME=`echo ${ARR[$i]} | awk -F "_" '{print $1}'`
printf "\n"
echo "[merging] $SAMPLE_NAME"
echo "samtools merge ./1-merged/$SAMPLE_NAME.bam ${ARR[$i]} ${ARR[$i+1]}"
samtools merge ./1-merged/$SAMPLE_NAME.bam ${ARR[$i]} ${ARR[$i+1]}
}
cd 1-merged
TOTAL_FILES=`find -iname '*.bam' | wc -l`
ARR=($(ls *.bam))
for ((i=0; i<$TOTAL_FILES; i+=1))
{
SAMPLE_NAME=`echo ${ARR[$i]} | awk -F "_" '{print $1}'`
printf "\n"
echo "[sorting] $SAMPLE_NAME "
echo "samtools view -b ${ARR[$i]} | samtools sort -o ./../2-sorted/$SAMPLE_NAME"
samtools view -b ${ARR[$i]} | samtools sort -o ./../2-sorted/$SAMPLE_NAME
}
cd ..
3.	Bring singlets back to merged ones.
4.	Calculate coverage after merging exactly as in step 6.3.
6.5	iVAR processing
Place in one directory
•	final (sorted and merged) BAM files
o	The following commands process pair of files at once. Thus, it should be two technical replicates per sample (two BAM files) only. No more, no less.
o	Rename all these replicates to have a similar beginning of the filename ending with “_“ symbol. E.g for 10A sample keep “10A_” as first symbols: 10A_r1.BAM and 10A_r2.BAM. 10A_.BAM and 10A_1.BAM are also fine. Any file like 10A._BAM or 10A.BAM will be wrong.
•	FASTA file (zika_primers.fa) with ZIKV primers which were used for targeted amplification of pool #1 and pool #2
•	FASTA file with the ZIKV reference genome (zikv-pr.fasta). The same one as in step 6.2.
•	TSV file with information on ZIKV primer pairs (zika_primer_pair_information.tsv). This file should be created manually in the text editor.
o	Example primer pair information file from iVAR manual (https://github.com/andersen-lab/ivar/blob/master/docs/MANUAL.md)
400_1_out_L    400_1_out_R
400_2_out_L    400_2_out_R
400_3_out_L    400_3_out_R
...
•	GFF file with information on ORFs (zikv-pr.GFF3).
o	To make GFF file, go to Genbank and find your reference genome (KU501215.1 for ZIKV).
o	Export GFF file: Send to: > Complete Record > Choose Destination: File > Format: GFF3 > Create file
o	Don not pay attention to the message “GFF file is not in GFF3 file format!” in the final iVar output. According to iVar author (https://github.com/andersen-lab/ivar/issues/23):
	The GFF3 file from NCBI has 4 lines starting with #! and those lines are showing the error. It's a simple fix and this should be done for next release. Irrespective of that though, ivar should proceed with translation as expected.
Note: Be very cautious because GFF works with ZIKV but did not work with PCV2. In PCV2 genome ORF2 transcribed in reverse direction and GFF gives misleading positions for AA.  Be careful working with viruses with genomes different from flaviviruses; for correct identification of AA positions use manual method. 
Run script that executes the following steps:
1.	Creating the main BAM file with primers. Creating the main BED file with primers.
2.	Counting and processing BAM files with replicates of sample reads one-by-one. BAM file with primers is not counted and is not processed in this pipeline.
e.	Making index.
f.	Trimming (soft clipping) primers from BAM.
g.	Sorting new BAM. Indexing BAM. Merging replicates to a merged BAM file.
h.	Processing merged BAM file representing both replicates.
i.	Making consensus FASTA file from merged BAM. Indexing it. Aligning primers to it and creating individual sorted BAM files with primers for each sample.
ii.	Creating the TSV file with masked primers. Creating individual BED file with primers.
i.	Removing reads with mismatches in primer sequences from BAM files representing replicates. Sorting of resulting BAM files.
j.	Generating the final TSV file with the list of mutations in each replicate. Only mutations above the threshold (3%; highlighted with red) are recorded.
k.	Filtering similar SNVs across replicates. Generating the final TSV file with the list of mutations in a sample (in each of both replicates).
•	Input files are highlighted with yellow. If you have different ones then edit the highlighted fragments accordingly.
clear
date
REFERENCE=zikv-pr.fasta
PRIMERS=zika_primers.fa
PRIMER_PAIRS=zika_primer_pair_information.tsv
GFF_FILE=zikv-pr.gff3
# We are counting original BAM files here because one extra BAM file will be created soon for primers.
TOTAL_FILES=`find -iname '*.bam' | wc -l`
ARR=($(ls *.bam))
bwa index $REFERENCE
bwa mem -k 5 -T 16 $REFERENCE $PRIMERS | samtools view -b -F 4 > primers.bam
bedtools bamtobed -i primers.bam > primers.bed
printf "\n\n"     
echo "Part 1: processing individual files with indexing, trimming, sorting, and indexing"
printf "\n"
for ((i=0; i<$TOTAL_FILES; i+=1)) {
    FILE_NAME=`echo ${ARR[$i]} | awk -F "." '{print $1}'`
    printf "\n"
    echo "[processing] $FILE_NAME"
    echo "[indexing] samtools index $FILE_NAME.bam"
    samtools index $FILE_NAME.bam
    echo "[trimming] ivar trim -b primers.bed -i ${ARR[$i]} -p $FILE_NAME.trimmed"
    ivar trim -b primers.bed -i ${ARR[$i]} -p $FILE_NAME.trimmed
    echo "[sorting] samtools sort $FILE_NAME.trimmed.bam -o $FILE_NAME.trimmed.sorted.bam"
    samtools sort $FILE_NAME.trimmed.bam -o $FILE_NAME.trimmed.sorted.bam
    echo "[indexing] samtools index $FILE_NAME.trimmed.sorted.bam"
    samtools index $FILE_NAME.trimmed.sorted.bam
}
printf "\n\n"
echo "Part 2: processing two replicates with merging, making and indexing consensus, primers processing (aligning, framing, and finding mismatches)"
printf "\n"
for ((i=1; i<$TOTAL_FILES; i+=2)) {
    SAMPLE_NAME=`echo ${ARR[$i]} | awk -F "_" '{print $1}'`
    FILE_NAME1=`echo ${ARR[$i-1]} | awk -F "." '{print $1}'`
    FILE_NAME2=`echo ${ARR[$i]} | awk -F "." '{print $1}'`
    printf "\n"
    echo "[processing] $SAMPLE_NAME"
    echo "[merging replicates] samtools merge $SAMPLE_NAME.merged.bam $FILE_NAME1.trimmed.sorted.bam $FILE_NAME2.trimmed.sorted.bam"
    samtools merge $SAMPLE_NAME.merged.bam $FILE_NAME1.trimmed.sorted.bam $FILE_NAME2.trimmed.sorted.bam
    echo "[making consensus] samtools mpileup -A -d 1000000 -Q 0 $SAMPLE_NAME.merged.bam | ivar consensus -p $SAMPLE_NAME.consensus"
    samtools mpileup -A -d 1000000 -Q 0 $SAMPLE_NAME.merged.bam | ivar consensus -p $SAMPLE_NAME.consensus
    echo "[primers processing: indexing consensus] bwa index $SAMPLE_NAME.consensus.fa"
    bwa index $SAMPLE_NAME.consensus.fa
    echo "[primers processing: conversion to sorted bam] bwa mem -k 5 -T 16 $SAMPLE_NAME.consensus.fa $PRIMERS | samtools view -bS -F 4 | samtools sort -o $SAMPLE_NAME.primers_consensus.bam"
    bwa mem -k 5 -T 16 $SAMPLE_NAME.consensus.fa $PRIMERS | samtools view -bS -F 4 | samtools sort -o $SAMPLE_NAME.primers_consensus.bam
    echo "[primers processing: creating BED] bedtools bamtobed -i $SAMPLE_NAME.primers_consensus.bam > $SAMPLE_NAME.primers_consensus.bed"
    bedtools bamtobed -i $SAMPLE_NAME.primers_consensus.bam > $SAMPLE_NAME.primers_consensus.bed
    echo "[primers processing: creating TSV file] samtools mpileup -A -d 1000000 --reference $SAMPLE_NAME.consensus.fa -Q 0 $SAMPLE_NAME.primers_consensus.bam | ivar variants -p $SAMPLE_NAME.primers_consensus -t 0.03"
    samtools mpileup -A -d 1000000 --reference $SAMPLE_NAME.consensus.fa -Q 0 $SAMPLE_NAME.primers_consensus.bam | ivar variants -p $SAMPLE_NAME.primers_consensus -t 0.03
    echo "[primers processing: masking] ivar getmasked -i $SAMPLE_NAME.primers_consensus.tsv -b $SAMPLE_NAME.primers_consensus.bed -f $PRIMER_PAIRS -p $SAMPLE_NAME.mismatches"
    ivar getmasked -i $SAMPLE_NAME.primers_consensus.tsv -b $SAMPLE_NAME.primers_consensus.bed -f $PRIMER_PAIRS -p $SAMPLE_NAME.mismatches
}
printf "\n\n"
echo "Part 3: processing individual files with removing reads with mismatches, sorting BAM-file, and calling SNVs"
printf "\n"
for ((i=0; i<$TOTAL_FILES; i+=1)) {
    SAMPLE_NAME=`echo ${ARR[$i]} | awk -F "_" '{print $1}'`
    FILE_NAME=`echo ${ARR[$i]} | awk -F "." '{print $1}'`
    printf "\n"
    echo "[processing] $FILE_NAME"
    echo "[removing reads] ivar removereads -i $FILE_NAME.trimmed.sorted.bam -p $FILE_NAME.trimmed.sorted.masked.bam -t $SAMPLE_NAME.mismatches.txt -b primers.bed"
    ivar removereads -i $FILE_NAME.trimmed.sorted.bam -p $FILE_NAME.trimmed.sorted.masked.bam -t $SAMPLE_NAME.mismatches.txt -b primers.bed
    echo "[sorting resulting BAM] samtools sort $FILE_NAME.trimmed.sorted.masked.bam -o $FILE_NAME.trimmed.sorted.masked.sorted.bam"
    samtools sort $FILE_NAME.trimmed.sorted.masked.bam -o $FILE_NAME.trimmed.sorted.masked.sorted.bam
    echo "[generating final TSV file] samtools mpileup -A -d 1000000 --reference $REFERENCE –B -Q 0 $FILE_NAME.trimmed.sorted.masked.sorted.bam | ivar variants -p $FILE_NAME.SNV -t 0.03 -r $REFERENCE -g GFF_FILE"
    samtools mpileup -A -d 1000000 --reference $REFERENCE –B -Q 0 $FILE_NAME.trimmed.sorted.masked.sorted.bam | ivar variants -p $FILE_NAME.SNV -t 0.03 -r $REFERENCE -g GFF_FILE
}
printf "\n\n"
echo "Part 4: Filtering the same SNVs from replicates"
printf "\n"
for ((i=1; i<$TOTAL_FILES; i+=2)) {
    SAMPLE_NAME=`echo ${ARR[$i]} | awk -F "_" '{print $1}'`
    FILE_NAME1=`echo ${ARR[$i-1]} | awk -F "." '{print $1}'`
    FILE_NAME2=`echo ${ARR[$i]} | awk -F "." '{print $1}'`
    printf "\n"
    echo "[processing] $SAMPLE_NAME"
    echo "[filtering SNVs across replicates] ivar filtervariants -p $SAMPLE_NAME.SNV_final $FILE_NAME1.SNV.tsv $FILE_NAME2.SNV.tsv"
    ivar filtervariants -p $SAMPLE_NAME.SNV_final $FILE_NAME1.SNV.tsv $FILE_NAME2.SNV.tsv
}
echo "[Sorting files to directories and renaming]"
mkdir 0.input && mv $PRIMERS $REFERENCE $PRIMER_PAIRS $GFF_FILE $_
mkdir 1.index_files && mv *.pac $_ && mv *.sa $_ && mv *.ann $_ && mv *.amb $_ && mv *.fai $_ && mv *.bai $_ && mv *.bwt $_
mkdir 1.Primers_info && mv primers.bam primers.bed $_
mkdir 2a.Consensus && mv *.consensus.fa $_
mkdir 2c.BED && mv *.bed $_
mkdir 2d.Consensus-TSV && mv *.primers_consensus.tsv $_
mkdir 3.Mismatches  && mv *.mismatches.txt $_
mkdir 4.Quality && mv *.qual.txt $_
mkdir 7.SNV-final && mv *.SNV_final.tsv $_
mkdir 6.SNV-singlets && mv *.SNV.tsv $_
mkdir 5.BAM-final && mv *.masked.sorted.bam $_
rm *.masked.bam
mkdir 2b.BAM-primers && mv *.primers_consensus.bam $_
rm *.merged.bam
rm *.sorted.bam
rm *.trimmed.bam
mkdir 0.BAM-original && mv *.bam $_
cd ./7.SNV-final
for FILE in *; do mv "$FILE" "${FILE/SNV_final./}"; done
cd ../5.BAM-final
for FILE in *; do mv "$FILE" "${FILE/trimmed.sorted.masked.sorted./}"; done
cd ..
date
echo "[FINISHED]"
6.6	Analysis of SNVs present in the positive control sample
Re-run the code provided in chapter 6.6 code for positive control samples. Change mutation detection threshold from 3% (0.03; highlighted with red color) to 0.01% (0.0001).
Check the resulting file with the threshold of 0.1% for the mutations detected in your samples with the threshold of 3%. Use reference values for your equipment to assess importance of the detected mutations in the positive control stock/inoculum.
For illumina Miseq, company states for 0.1% (doi: 10.1111/j.1755-0998.2011.03024.x), while researchers report 0.24% (doi: 10.1038/s41598-018-29325-6) and 0.30-0.46% (doi: 10.1186/gb-2013-14-5-r51) as an error rate.
6.7	Calculate coverage for the final masked and sorted BAM files
Calculate coverage in final BAM files generated by iVAR exactly as in step 6.3.
6.8	SNVs processing
Import final TSV files to Microsoft Excel for Microsoft Windows (this will not work in Libre Office).
1.	Put your TSV files in a dedicated folder.
2.	Open Excel. Create new empty document.
3.	Data > New Query > From File > From Folder > <Provide the pathway for your folder with TSV files> OK > Combine > Combine&Edit > Change “Data Type Detection” to “Do not detect data types” and “Delimiter” to “Tab”> OK > File > Close&Load
a.	At this step, mutations with any coverage are imported.
4.	For each technical replicate delete manually lines representing insertions/deletions (i.e., mutations labeled with plus or minus sign. E.g. column E contains “+A”).
5.	For each technical replicate delete manually lines representing mutations with coverage <400 (Columns R and AB: TOTAL_DP).
6.	Add the “Group” column and fill it with information on groups.
Optional. Use Microsoft Excel to find the masked regions (coverage <400) in the positive control sample and exclude mutations in all samples falling in these regions. Information about masked regions is available in the report of the R script generated in step 6.6.
