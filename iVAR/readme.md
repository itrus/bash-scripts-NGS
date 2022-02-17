# Single nucleotide variation analysis in Zika virus (ZIKV) and porcine circovirus 2 (PCV2) genomes

The current manual was published as Supplementary File in Viruses. But as soon as the code can not be updated there, I decided to upload all the scripts to GitHub.

> Udenze, Daniel, Ivan Trus, Henry Munyanduki, Nathalie Berube, and Uladzimir Karniychuk. 2021. "The Isolated in Utero Environment Is Conducive to the Emergence of RNA and DNA Virus Variants" _Viruses_ 13, no. 9: 1827. https://doi.org/10.3390/v13091827

## 0 SUMMARY

This manual describes a bioinformatics pipeline used to identify SNVs (single nucleotide variants) in Next-Generation Sequencing (NGS) data from ZIKV and PCV2 genomes. The pipeline can be adapted to analyse NGS data from other viruses. In the first step, paired-end reads from Illumina NGS data are preprocessed. After confirming the expected coverage and depth in the PhiX sequencing control, SNVs are detected with iVAR package.


![overview pre-processing](https://user-images.githubusercontent.com/9166776/153178845-8a7f503e-77e5-4d4b-a313-1b030ea1ea11.png)
Figure 1. Preprocessing

![overview variant calling](https://user-images.githubusercontent.com/9166776/153178924-816a78f3-cd8d-49e1-b730-7c60ac2edf11.png)
Figure 2. SNVs calling with iVAR

## 1 PREREQUISITES

We use Ubuntu Linux for data preprocessing. Run all commands in the console/terminal/shell. Configure your Linux shell. Go to the Terminal preferences > Unnamed profile > Scrolling > Remove the checkmark from “Limit scrollback to ...”.
For running each command, you need to enter the proper directory with the terminal session. For example, for trimming PhiX, you will need to use the directory /home/ngs1/zika virus in utero heterogeneity IP versus IC/temporary. To enter it with terminal, execute the following command:
cd "/home/ngs1/zika virus in utero heterogeneity IP versus IC/temporary"

After running each command:

- If several commands were copy-pasted, check that the last command was executed. If not, press Enter to launch the last command.
-	Read all output in Terminal to check if there are no errors.
-	Copy-paste the terminal session's output to a file (log.txt) and save it in the working folder for each step.

For each step, you need an empty working directory. Bring input files; after successful execution, delete input files, but preserve all newly created files (results) and reference files (e.g., GENOME-, BED-files, etc.). This allows us to go back and to repeat any step under the same or different conditions, if necessary.


## 2 PREPARE THE LINUX ENVIRONMENT

The step-by-step installation manual for iVAR is available online (https://andersen-lab.github.io/ivar/html/installpage.html).

First, install all required dependencies:
```
sudo apt autoremove
sudo apt install autotools-dev gcc zlib1g-dev libbz2-dev \
liblzma-dev make libncurses5-dev autoconf g++ bwa r-base \
bedtools parallel default-jre mc krusader krename kate gawk
```
To install HTSlib libraries, download HTSlib 1.14 (http://www.htslib.org/download/), then finish the installation procedure.
```
tar -xvf ./htslib-1.14.tar.bz2
cd ./htslib-1.14/
./configure
make
sudo make install
```
If HTSlib is installed in a non-standard location (e.g., **/usr/local/lib/**) and iVAR cannot find it, please add the following line to your BASH configuration file (**~/.bash_profile or ~/.bashrc**) so that iVar can find HTSlib dynamic libraries during runtime. Restart is not needed.

```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/
```

To install SAMlib libraries, download SAMtools 1.14 (http://www.htslib.org/download/) then finish the installation procedure.
```
tar -xvf ./samtools-1.14.tar.bz2
cd ./samtools-1.14/
./configure
make
sudo make install
```

Install iVAR (https://github.com/andersen-lab/ivar) with running the **install-ivar.sh** script.

Install Trimmomatic (v 0.40; binary; http://www.usadellab.org/cms/?page=trimmomatic). Download and unpack it to /home/user/bin/Trimmomatic-0.40/. UGENE on Windows contains a build-in version of Trimmomatic and can replace the Linux version.
```
unzip ./Trimmomatic-0.40.zip
```
To install RStudio, download from the official website (https://rstudio.com/products/rstudio/download/#download) version for Ubuntu 18 (or your version). Then install the downloaded file (Right mouse button > Open With Software Install > Install)

## 3 QUALITY CONTROL

### 3.1 Illumina Seq Analysis Viewer

Use Illumina Seq Analysis Viewer 2.4.7.0 for Windows to see the general report of the sequencing (Browse > Use MiSeq output whole folder as input for Illumina Seq Analysis Viewer). The full manual for the software package is available on Illumina’s website (https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/sav/sequencing-analysis-viewer-user-guide-15020619-f.pdf).

Create a report in Word containing charts and data. To enlarge the chart click the chevron in the upper right corner of each chart. To copy-paste charts, click it with the right mouse button and select Copy To Clipboard.
1. Analysis tab
   - Data by Cycle
     - Error Rate (percentage of wrong bases in the PhiX control reads).
     - Called Int  (corrected intensity of base calls)
   - QScore Distribution (Q20 corresponds to 1% error rate, Q30 corresponds to 0.1% error rate)
   - Data by Lane
     - % Aligned (Frequency of reads representing PhiX)
   - Qscore Heatmap (Q-scores by cycle)
2. Summary tab. Press Copy To Clipboard to copy-paste data.
   - Total – Yield Total (G) (Total number of Gbp saved)
   - Total – Aligned (%) (Frequency of reads representing PhiX)
   - Error Rate (%) (percentage of wrong bases in the PhiX control reads)
3. Indexing tab.
   - Total reads (total number of reads saved)
   - % reads Identified (PF) (percentage of reads with barcodes)
   - Min and Max (percentage of the least and most frequent barcodes)
   - Click on a table with individual indexes with the left mouse button, then press Control+A > Control+C to copy-paste the table content.
   - Barchart contains a graphical representation of the table. Copy-paste chart by clicking it with the right mouse button and selecting Copy To Clipboard.

Use the corresponding folder to save the results (screenshots).

### 3.2 FastQC

Alternatively, you can download and process FASTQ.GZ files with FastQC.
1. Download FastQC (https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc) and unpack it to /home/user/bin.
2. Run qc.sh script.

## 4 Trimmomatic

Copy-paste to a new working folder FASTQ.GZ files containing paired reads and trimmomatic.sh script. Launch the script to proceed with Trimmomatic.
- This workflow is only valid for the paired-ended reads (PE). You cannot use it for single-ended (SE) reads.
- Default workflow parameters: ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36. Adapt it in the trimmomatic.sh script if necessary.
- Illimina HiSeq is using Phred33 quality scores and needs TruSeq3-PE-2.fa file from the Trimmomatic package for removing adaptors. Input files for adaptors may vary depending on the machine, etc. Check it each time.
- Phred >20 parameter used for quality trim corresponds to 99% accuracy in nucleotides.
- Adapt the directory pathway for the Trimmomatic package in the trimmomatic.sh script if necessary.
- After trimming, check the report of Trimmomatic carefully for errors and save it as a log file (log.txt). The output is 4 files sorted into two directories:
  - 2 files with paired reads (forward and reverse): <SAMPLE>_1P.fastq.gz and <SAMPLE>_2P.fastq.gz stored in ./trimmed/paired directory.
  - 2 files with unpaired reads (forward and reverse): <SAMPLE>_1U.fastq.gz and <SAMPLE>_2U.fastq.gz stored in ./trimmed/unpaired directory.
- Collect and keep generated output files from the newly created folder named ‘trimmed’. For the following steps, you will need only 2 files with paired reads.
   
## 5 Align reads to the reference with BWA

- Make a copy of output FASTQ.GZ files from Trimmomatic corresponding to paired reads (2 files per sample). Don’t mix paired and unpaired files after Trimmomatic. Proceed only with paired ones. For the PhiX control sample name could be "Undetermined".
Provide a file representing the reference for the alignment.
- This is an important step! Get the CORRECT reference (FASTA file). Place the reference FASTA file in the same folder as FASTQ.GZ files from Trimmomatic. Check the bwa.sh script and edit it if necessary. You need to have the correct reference mentioned/uncommented in the first lines of the script.
  - For example, for PhiX control use the Illumina website (https://support.illumina.com/sequencing/sequencing_software/igenome.html). Genbank has the PhiX sequence [Genbank: NC_001422.1 that is 5 nt different from the Illumina’s one.
   - For the PR strain of ZIKV I downloaded FASTA file (ZIKV-PR.fasta) from the Genbank [Genbank: KU501215.1].
   - For the BaPCV2b strain of PCV2 I downloaded FASTA file (BaPCV2b.fasta) from the Genbank [Genbank: FJ233905.1].
- Run bwa.sh script to align paired reads to the reference. First index files for the reference genome will be created. Then as the resulting output sorted and merged BAM file named as <SAMPLE>.BAM will be created. Unmapped reads are skipped (“F” parameter). After successfully run, delete the copy-pasted input FASTQ.GZ files to save the space and other tiny files (e.g. index), but keep the generated output (BAM file) and logging information from the shell/terminal.

## 6 Coverage analysis
### 6.1 UGENE
To visually inspect the coverage, open the resulting BAM file with UGENE.
Control (e.g. PhiX genome) should be covered fully (100% genome coverage) with a high coverage depth (>10’000x).

### 6.2 Script

For batch processing of all samples at once and detailed analysis:
1. Create a GENOME file from the reference FASTA file (e.g., ZIKV-PR.fasta). 
```
INDEX=ZIKV-PR.fasta
samtools faidx $INDEX
awk {'printf $1 "\t" $2'} $INDEX.fai > $INDEX.genome
```
2. Copy the final (merged and sorted) BAM files to a new directory. Add to the same folder the correct GENOME file (e.g. ZIKV-PR.fasta.genome). The following commands will produce TXT files with the coverage. Besides this BED files are created. You may delete them as soon as they will not be required further.
```
find ./ -iname '*.bam' | parallel -j 4 --progress 'bedtools bamtobed -i {} > {.}.bed | echo {/}'
find ./ -iname '*.bed' | parallel -j 4 --progress 'bedtools genomecov -d -i {} -g ZIKV-PR.fasta.genome > {.}.txt'
```
3. Place the R script (coverage.Rmd) and TXT files with coverages in a separate (new) directory. Make sure that only the relevant TXT files were placed in this folder. Analyse resulting files with the coverage.Rmd R script.
   - Open and run the R script with RStudio.
   - Use **Code > Run Region > Run All (Control+Alt+R)** if you want to see results in RStudio.
   - Use **File > Knit Document** if you need a report (HTML file) generated. This file can be opened with any web browser. To create a report in PDF file format, you need first to install the LaTeX distribution. 

## 7 Merging BAM files representing one sample but pool #1 and pool #2
Do this step when DNA representing the whole ZIKV or PCV2 genome was amplified with two different primer pools (see materials and methods for PrimalSeq NGS) and was not mixed during NGS library preparation. If pools were mixed, each biological sample is represented by only one BAM file and merging is not required. Pools are not the same as technical replicates.
1. Rename all files manually. E.g. for Pool #1 and Pool #2: files 4A.bam and 5A.bam should be renamed to 4A.bam and 4A_5A.bam.
2. The script mergeBAMs.sh processes a pair of files at once. Thus, temporarily remove singlets (negative control samples represented by pool #1 or pool #2 only) to a separate folder.
3. Calculate coverage after merging exactly as in step 6.

## 8 iVAR processing
Place in one directory:
1. Final BAM files
   - The following commands process pair of files at once. Thus, it should be two technical replicates per sample (two BAM files) only. No more, no less. Rename all these replicates to have a similar beginning of the filename ending with the “_“ symbol. E.g for the 10A sample keep “10A_” as the first symbols: 10A_r1.BAM and 10A_r2.BAM. 10A_.BAM and 10A_1.BAM are also fine. Any file like 10A._BAM or 10A.BAM will be wrong.
2. FASTA file (primers.fasta) with ZIKV primers which were used for targeted amplification of pool #1 and pool #2
3. FASTA file with the ZIKV reference genome (ZIKV-PR.fasta). The same one as in steps 5 and 6.
4. TSV file with information on ZIKV primer pairs (primer_pair_information.tsv). This file should be created manually in the text editor.
   - Example primer pair information file from iVAR manual (https://github.com/andersen-lab/ivar/blob/master/docs/MANUAL.md)
```
400_1_out_L    400_1_out_R
400_2_out_L    400_2_out_R
400_3_out_L    400_3_out_R
```
5. GFF file with information on ORFs (ZIKV-PR.gff3) for AA translations.
   - To make GFF file, go to Genbank and find your reference genome (KU501215.1 for ZIKV). Export GFF file: **Send to: > Complete Record > Choose Destination: File > Format: GFF3 > Create file**
   - Do not pay attention to the message "_GFF file is not in GFF3 file format!_" in the final iVar output. According to iVar author (https://github.com/andersen-lab/ivar/issues/23): _The GFF3 file from NCBI has 4 lines starting with #! and those lines are showing the error. It's a simple fix and this should be done for next release. Irrespective of that though, ivar should proceed with translation as expected._
   - Note: Be very cautious because GFF works with ZIKV but did not work with PCV2. In PCV2 genome ORF2 is transcribed in reverse direction and GFF gives misleading positions for AA. Be careful working with viruses with genomes different from flaviviruses; for correct identification of AA positions use the manual method. 
6. run-ivar.sh script
   - Input files are mentioned in the header of the file. Check the filenames carefully. If you have different ones then edit the script accordingly.
   - The THRESHOLD parameter in the header of the file has a default setting of 3% (0.03). It stands for the minimal level of SNV prevalence in both replicates that will be reported. You may re-run iVAR for positive control or all samples after decreasing the mutation detection threshold from 3% (0.03) to 0.01% (0.0001). Check the resulting file representing the control (virus stock or inoculum) with the threshold of 0.1% for the mutations detected in your samples with the threshold of 3%. You may use reference values for your equipment to assess the importance of the detected mutations in the positive control stock/inoculum. For Illumina Miseq, company states for 0.1% (doi: 10.1111/j.1755-0998.2011.03024.x), while researchers report 0.24% (doi: 10.1038/s41598-018-29325-6) and 0.30-0.46% (doi: 10.1186/gb-2013-14-5-r51) as an error rate.

Execute the run-ivar.sh script that goes through the following steps:
1. Creating the main BAM file with primers. Creating the main BED file with primers.
2. Counting and processing BAM files with replicates of sample reads one by one. BAM file with primers is not counted and is not processed in this pipeline.
   - Making index.
   - Trimming (soft clipping) primers from BAM.
   - Sorting new BAM. Indexing BAM. Merging replicates to a merged BAM file.
   - Processing merged BAM file representing both replicates.
     - Making consensus FASTA file from merged BAM. Indexing it. Aligning primers to it and creating individual sorted BAM files with primers for each sample.
     - Creating the TSV file with masked primers. Creating individual BED files with primers.
   - Removing reads with mismatches in primer sequences from BAM files representing replicates. Sorting of resulting BAM files.
   - Generating the final TSV file with the list of mutations in each replicate. Only mutations above the threshold (3%; highlighted with red) are recorded.
   - Filtering similar SNVs across replicates. Generating the final TSV file with the list of mutations in a sample (in each of both replicates).

## 9 Final steps
1. Calculate coverage in final BAM files generated by iVAR exactly as in step 7.
2. Combine all SNVs in a spreadsheet
   - To import final TSV files to Microsoft Excel for Microsoft Windows (this will not work in the same way in Libre Office) put your TSV files in a dedicated folder.
   - Open Excel. Create a new empty document.
   - **Data > New Query > From File > From Folder > <Provide the pathway for your folder with TSV files> OK > Combine > Combine&Edit > Change “Data Type Detection” to “Do not detect data types” and “Delimiter” to “Tab”> OK > File > Close&Load**. At this step, mutations with any coverage are imported to Excel from all provided TSV files at once.
   - For each technical replicate delete manually lines representing insertions/deletions (i.e., mutations labeled with plus or minus sign. E.g. column E contains “+A”).
   - For each technical replicate delete manually lines representing mutations with coverage <400 (Columns R and AB: TOTAL_DP).
   - Add the “Group” column and fill it with information on groups.
   - Optional step. Use Microsoft Excel to find the masked regions (coverage <400) in the positive control sample and exclude mutations in all samples falling into these regions. Information about masked regions is available in the report of the R script generated in step 7.
