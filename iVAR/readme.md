# Single nucleotide variation analysis in Zika virus (ZIKV) and porcine circovirus 2 (PCV2) genomes

## 1	SUMMARY

This manual describes a bioinformatics pipeline used to identify SNVs (single nucleotide variants) in Next-Generation Sequencing (NGS) data from ZIKV and PCV2 genomes. The pipeline can be adapted to analyse NGS data from other viruses. In the first step, paired-end reads from illumina NGS data are preprocessed. After confirming the expected coverage and depth in the PhiX sequencing control, SNVs are detected with iVAR package.
The current manual was published as Supplementary File in Viruses. But as soon as code can not be updated I decided to upload all the scripts to GitHub.

> Udenze, Daniel, Ivan Trus, Henry Munyanduki, Nathalie Berube, and Uladzimir Karniychuk. 2021. "The Isolated in Utero Environment Is Conducive to the Emergence of RNA and DNA Virus Variants" _Viruses_ 13, no. 9: 1827. https://doi.org/10.3390/v13091827

## 2	PREREQUISITES

We use Ubuntu Linux for data preprocessing. Run all commands in the console/terminal/shell. Configure your Linux shell. Go to the Terminal preferences > Unnamed profile > Scrolling > Remove the checkmark from “Limit scrollback to ...”.
For running each command, you need to enter the proper directory with the terminal session. For example, for trimming PhiX, you will need to use the directory /home/ngs1/zika virus in utero heterogeneity IP versus IC/temporary. To enter it with terminal, execute the following command:
cd "/home/ngs1/zika virus in utero heterogeneity IP versus IC/temporary"

After running each command:

o	If several commands were copy-pasted, check that the last command was executed. If not, press Enter to launch the last command.

o	Read all output in Terminal to check if there are no errors.

o	Copy-paste the terminal session's output to a file (log.txt) and save it in the working folder for each step.

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

Install iVAR (https://github.com/andersen-lab/ivar).

```
wget https://github.com/andersen-lab/ivar/archive/refs/heads/master.zip
unzip ./master.zip
cd ./ivar-master/
./autogen.sh
./configure
make
sudo make install
```

Test installation of iVAR 1.3.1.

```
ivar version
```

Install Trimmomatic (v 0.40; binary; http://www.usadellab.org/cms/?page=trimmomatic). Download and unpack it to /home/user/bin/Trimmomatic-0.40/. UGENE on Windows contains a build-in version of Trimmomatic and can replace the Linux version.
```
unzip ./Trimmomatic-0.40.zip
```
To install RStudio, download from the official website (https://rstudio.com/products/rstudio/download/#download) version for Ubuntu 18 (or your version). Then install the downloaded file (Right mouse button > Open With Software Install > Install)
