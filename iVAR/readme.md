# Single nucleotide variation analysis in Zika virus (ZIKV) and porcine circovirus 2 (PCV2) genomes

The current manual was published as Supplementary File in Viruses. But as soon as code can not be updated there, I decided to upload all the scripts to GitHub.

> Udenze, Daniel, Ivan Trus, Henry Munyanduki, Nathalie Berube, and Uladzimir Karniychuk. 2021. "The Isolated in Utero Environment Is Conducive to the Emergence of RNA and DNA Virus Variants" _Viruses_ 13, no. 9: 1827. https://doi.org/10.3390/v13091827

## 1	SUMMARY

This manual describes a bioinformatics pipeline used to identify SNVs (single nucleotide variants) in Next-Generation Sequencing (NGS) data from ZIKV and PCV2 genomes. The pipeline can be adapted to analyse NGS data from other viruses. In the first step, paired-end reads from illumina NGS data are preprocessed. After confirming the expected coverage and depth in the PhiX sequencing control, SNVs are detected with iVAR package.


![Figure 1. Preprocessing](https://github.com/itrus/bash-scripts-NGS/raw/main/iVAR/overview%20pre-processing.png)
Figure 1. Preprocessing

![Figure 2. Variant calling](https://raw.githubusercontent.com/itrus/bash-scripts-NGS/main/iVAR/overview%20variant%20calling.png)
Figure 2. SNVs calling with iVAR

## 2	PREREQUISITES

We use Ubuntu Linux for data preprocessing. Run all commands in the console/terminal/shell. Configure your Linux shell. Go to the Terminal preferences > Unnamed profile > Scrolling > Remove the checkmark from “Limit scrollback to ...”.
For running each command, you need to enter the proper directory with the terminal session. For example, for trimming PhiX, you will need to use the directory /home/ngs1/zika virus in utero heterogeneity IP versus IC/temporary. To enter it with terminal, execute the following command:
cd "/home/ngs1/zika virus in utero heterogeneity IP versus IC/temporary"

After running each command:

- If several commands were copy-pasted, check that the last command was executed. If not, press Enter to launch the last command.
-	Read all output in Terminal to check if there are no errors.
-	Copy-paste the terminal session's output to a file (log.txt) and save it in the working folder for each step.

For each step, you need an empty working directory. Bring input files; after successful execution, delete input files, but preserve all newly created files (results) and reference files (e.g., GENOME-, BED-files, etc.). This allows us to go back and to repeat any step under the same or different conditions, if necessary.


## 3	PREPARE THE LINUX ENVIRONMENT

The step-by-step installation manual for iVAR is available online (https://andersen-lab.github.io/ivar/html/installpage.html).

First install all require dependencies:
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

## 4	QUALITY CONTROL

Use Illumina Seq Analysis Viewer 2.4.7.0 for Windows to see the general report of the sequencing (Browse > Use MiSeq output whole folder as input for Illumina Seq Analysis Viewer). The full manual for the software package is available at illumina’s website (https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/sav/sequencing-analysis-viewer-user-guide-15020619-f.pdf).

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
   - Click on a table with individual indexes with left mouse button, then press Control+A > Control+C to copy-paste the table content.
   - Barchart contains a graphical representation of the table. Copy-paste chart by clicking it with the right mouse button and selecting Copy To Clipboard.

Use the corresponding folder to save the results (screenshots).
