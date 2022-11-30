# Bioinformatics pathway for RNA-Seq
## 1 SUMMARY
This manual describes the bioinformatics to identify DEGs (differentially expressed genes) and affected GOs (gene ontology pathways) in RNA-seq data. For training purposes, starting files are SRA files downloaded from the SRA database. In the first step they are converted to FASTQ.GZ files (Figure 1).

![rna-seq](https://user-images.githubusercontent.com/9166776/155319981-02d5fe89-bb44-4a98-b552-8a2cbfdfbd40.png)

Figure 1. Processing pathway.

If your starting data are paired-end reads from illumina NGS and contain adaptor sequences, please process your data with Trimmomatic to remove the adaptors (see Apendix). And then proceed to Kallisto analysis. If your starting data is made of single-end reads, the Kallisto Quantify transcripts (see below, 5.2) should be adjusted. Because Kallisto is expecting to get two FASTQ files for each sample/animal, while single-reads will give only one per sample.
In general, Kallisto is (a) making index for porcine genome (Sus scrofa, version 11.1) and (b) counting transcripts. In the next step, the R-script is downloading from BioMart the database of gene names corresponding to your transcripts . Then another R script is performing analysis of DEGs and GO pathways based on existing experimental factors/treatments/groups.
There are two R scripts in the pipeline. The first script has one parameter that may be edited: reference genome used for the whole analysis. In the second script, adjustments will depend on experimental factors/treatments taken into analysis. 

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
For each analysis step you need an empty working directory. Bring input files, then after successful processing delete input files, but leave all newly created files (results) and reference files (e.g. IDX-files, etc.). This allows to go back and to repeat any step starting with any initial point. Besides this, it comes possible to investigate differences between files created under different parameters/conditions.
- Create subfolders required for storing results of each procedure.

## 3 PREPARE THE LINUX ENVIRONMENT
First install all the require dependencies:
```
sudo apt autoremove
sudo apt install parallel gawk gzip default-jre
```
During installation, machine will ask admin's password.
Depending from machine and changes in packages, it is possible that other installations, e.g., other packages and new version of R, will be required.
Copy-paste commands from Word to Konsole doesn't always work (e.g., if there are invisible symbols from Word). If Konsole gives strange message – e.g., tar invalid option –“^” – then retype the same comman manually.

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

### SRA toolkit
Install the SRA Toolkit to work with SRA files (This step is necessary if RNA-seq data from the SRA database are used, e.g., for training or re-analysis of published dat. In a new experiment, raw files will come directly from the illumina machine). To install the SRA Toolkit:
1. Go to https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software and download	the file called “Ubuntu Linux 64 bit architecture”.
2. Unpack it
```
tar –xzf sratoolkit.2.10.8-ubuntu64.tar.gz
```
3. Move a directory with unpacked sratoolkit into ~/bin (symbol “~” stands for your home directory. E.g. /home/user/)
4. Check if it works well.
```
~/bin/sratoolkit.2.10.8-ubuntu64/bin/prefetch -V
```

### Kallisto
Install Kallisto:
1. Download it from http://pachterlab.github.io/kallisto/download 
2. Unpack it
```
tar -xvf kallisto_linux-v0.46.1.tar.gz
```
3. Move it to ~/bin
4. Test it
```
~/bin/kallisto/kallisto
```

### Cytoscape
Install Cytoscape with the required plugins:
1. Download Cytoscape (https://cytoscape.org/download.html).
2. Run commands to lunch the installer.
```
chmod +x Cytoscape_3_8_1_unix.sh
./Cytoscape_3_8_1_unix.sh
```
4. Follow the dialog: accept everything and select the pathway for the program.
5. Launch it using the icon on the Desktop.
6. Install several plug-ins.
   - **Apps > App manager**
   - Install 4 plug-ins: “EnrichmentMap”, “clusterMaker2”, “WordCloud”, “AutoAnnotate” or 4-in-1 plug-in “EnrichmentMap Pipeline Collection”.

Cytoscape analysis can be done in Windows. Use the same link (in Windows it will automatically download the Windows version) and follow the same instructions.

## 4 MAKING FASTQ.GZ, ACCESSION LIST, AND METADATA FILES
### 4.1 Download the Accession list – this section is for training/re-analyzing of publicly available data (see 4.3 for new data)
To make the Accession list with all of the samples visit the SRA Run Selector (https://trace.ncbi.nlm.nih.gov/Traces/study/?) and use the correct accession number for your study to find the page in the database with all the samples. E.g. use:
- PRJNA407675 for Pigs model for ZIKV (_Darbellay, J., Cox, B., Lai, K., Delgado-Ortega, M., Wheler, C., Wilson, D., ... & Karniychuk, U. (2017). Zika virus causes persistent infection in porcine conceptuses and may impair health in offspring. EBioMedicine, 25, 73-86._)
- PRJNA573521 for Subclinical ZIKV in pigs (_Trus, Ivan, et al. "Subclinical in utero Zika virus infection is associated with interferon alpha sequelae and sex-specific molecular brain pathology in asymptomatic porcine offspring." PLoS pathogens 15.11 (2019): e1008038._)
Go to **Select** tab. Then find the line with all samples (**Total**) and column with download links (**Download**). Download **Accession List** with all the samples (**SraAccList.txt**). Download also Metadata with data on the experiment setup and factors (**SraRunTable.txt**). This file will be necessary at step 6.1 for sample, group and factors identification.
### 4.2 Download and process SRA files – this section is for training/re-analyzing of publicly available data (see 4.3 for new data)
This protocol is for paired reads only. As soon as we have paired reads, two FASTQ.GZ files are generated from each SRA file.
This process requires a lot of processor time and of free disk space. The current computer (Kubuntu 20.04 on Intel Core i5-4590 with 32 Gb RAM) converts **SRA > FASTQ > FASTQ.GZ** with the speed of 5.9 Gb/h. This conversion process requires a lot of storage space: 7.2x-8.4x. E.g. to covert 10 Gb of SRA files you need to provide 84 Gb of empty space. After the end of the process you can delete the original SRA files to save the disk space.
Following commands using the **Accession list** will (a) download and (b) validate/verify corresponding SRA files. Then it will (c) convert them to FASTQ file-format and (d) compress to GZ file-format. Adapt directory pathways for the SRA toolkit and filename for the **Accession file** if necessary.
```
date
cat SraAccList.txt | xargs -i ~/bin/sratoolkit.2.10.8-ubuntu64/bin/prefetch -O ./ {}
date
find ./ -iname '*.sra' | parallel -j 4 --progress '~/soft/sratoolkit.2.10.8-ubuntu64/bin/vdb-validate {}'
date
find ./ -iname '*.sra' | parallel -j 4 --progress '~/soft/sratoolkit.2.10.8-ubuntu64/bin/fastq-dump --split-files {}'
date
find ./ -iname '*.fastq' | parallel -j 4 --progress 'gzip -v9 {}'
date
```
### 4.3 Create Accession list and Metadata files
If data obtained from a new illumina run, create Metadata (e.g./reference/RunTable.txt).
In Excel create columns and fill data in:
- SampleID - name of samples should match the existing illumina output fastq.gz files;
- Treatment – for example, ZIKV and control; 
And so on, create enough of columns to reflect all the experimental factors and experimental subgroups. 
Save this Excel file as .csv file; then, manually change extension to .txt.
This table will be used to compare gene expression in different experimental subgroups, see 6.2.  

## 5 KALLISTO
### 5.1 Build pig transcriptome index
In this experiment, we use the pig (_Sus scrofa_) reference genome. Download the latest reference genome for your experimental organism.

1. Download the reference pig transcriptome (cDNA and GTF files).
   - This is an important step! Get the correct reference. And starting from this step and until the end of this protocol work with only one reference genome at a time. There are plenty of different pig breeds sequenced and they are being updated regularly. Stay with only one reference genome within one experiment.
   - E.g. for getting porcine genome (Sscrofa11.1) visit http://uswest.ensembl.org/Sus_scrofa/Info/Index?db=core, click **Download FASTA**, click **cdna**, and download corresponding cDNA, e.g.,  _Sus_scrofa.Sscrofa11.1.cdna.all.fa.gz_ . Visit http://uswest.ensembl.org/Sus_scrofa/Info/Index?db=core, click **Download GTF**, and download corresponding GTF file, e.g.,  _Sus_scrofa.Sscrofa11.1.101.gtf.gz_ .
   - E.g. for the recent human genome, visit http://www.ensembl.org/Homo_sapiens/Info/Index, , click **Download FASTA**, click **cdna**, and download corresponding cDNA, e.g.,  _Homo_sapiens.GRCh38.cdna.all.fa.gz_ . Visit http://uswest.ensembl.org/Sus_scrofa/Info/Index?db=core, click **Download GTF**, and download corresponding GTF file, e.g.,  _Homo_sapiens.GRCh38.105.gtf.gz_ .
2. Make index for the genome (IDX file) from downloaded cDNA file and validate it with GTF file wih the **index.sh** script. This is a fast process. It takes ~5 minutes. Check the header of the script and modify the header lines with filenames if necessary. BED and IDX files will be created. You need only the IDX file to continue.

### 5.2 Quantify the transcripts
1. Bring to a new directory
   - All paired FASTQ.GZ or FASTQ files
   - Index file
   - GTF file
   - the **quant.sh** script. Modify filenames in the header of the file if necessary.
2. Execute the **quant.sh** script. The process goes relatively fast (11.6 Gb/h).
3. Resulting output is sorted in several directories (one directory for each sample). Do not rename directories or files. Do not move files to a single directory. Preserving the original file arrangement is important for proper import with **tximport** in the following analysis with the R script.
   - JSON file with run information
   - H5 file with binary data on transcript abundance in HDF5 file format
   - TSV file with plain-text information on transcript abundance

## 6 ANALYSIS IN R
Analysis in R is broken into two parts:
The first R script (RNA-seq-initial.Rmd) downloads gene names from the reference Ensembl genome, in this experiment from Sus Scrofa reference genome, and connects them to our experimental RNA-seq transcripts that were quantified by Kallisto at step 5.1, and generates the output file with all Gene expression data (**GE.Rdata**).
The second script (RNA-seq-final.Rmd) doesn’t need any connection to the Internet. It takes results of the first script (**GE.Rdata**) and analysis experimental RNA-seq data using the provided model that can be modified (e.g., comparing gene expression between Control and ZIKV-infected groups; include fetal sex as a factor into analysis, etc.).
### 6.1 The first R script is preparing data for analysis
This script (RNA-Seq-initial.Rmd) creates database of reference genes and links it to the experimental RNA-seq transcripts.
**This is an important step!** Get the correct reference: the same animal and the same version and the same breed (e.g. Sus scrofa 11.1). The database of genes and their respective transcripts should be created once (file named **Tx2Gen.Rdata**) and loaded later if necessary. Because it should be the same version as one used for Kallisto above (step 5.1).
The current script is searching for the most recent porcine genome (parameter _SearchString <- "sscrofa"_). Adjust it if necessary (e.g. for human genome it should be  _SearchString <- "GRCh38.p13"_). Check the output of the R script to see if it was working as expected.
Then RNA-Seq-initial.Rmd script loads all experimental Kallisto results and aligns experimental RNA-seq transcripts to the reference genome database. To make sure that identification of experimental RNA-seq transcripts and identification of genes in the reference dataset do match, you may go to the ENSEMBL website (http://uswest.ensembl.org/Sus_scrofa/Info/Index?db=core). Find 2-3 transcripts manually (copy-paste ID of the transcript), search and check if correct gene names were assigned properly.
Then script converts absolute counts/reads to relative counts-per-million reads (CPM) and cleans the dataset. For cleaning we keep Genes only with >1 CPM in at least 3 samples. Some people keep threshold at 0.2 CPM. Smarter approach is to adapt the threshold to each sample individually (https://doi.org/10.12688/f1000research.9005.3).
1. Place to a separate (new) directory
   - the first R script (RNA-Seq-initial.Rmd)
   - folder containing subfolders with transcript counts for each sample (“./kallisto”) (./ - means working directory; the actual folder name is – **kallisto**). It is important to have no extra files added by user to this folder; the workflow fails often in such cases.
   - file with Metadata for the whole trial (RunTable.txt)
2. Open and run the R script with R-Studio.
   - Use **Code > Run Region > Run All** (Control+Alt+R) if you want to see results in R-Studio.
   - Use **File > Knit Document** if you need a report file to be generated.
3. Script creates three new files:
   - RNA-Seq-initial.html – HTML report;
   - Tx2Gen.Rdata – the database of reference genes and link to their corresponding transcripts;
   - GE.Rdata contains all the data on Gene Expression and sample arrangement to groups (Factors; e.g, Control or ZIKV, etc.). This file will be an input for the next script.

### 6.2 The second R script makes the analysis
While the previous R script needs only one parameter to be adjusted (the reference genome), the second script has multiple parameters to consider.
1. Place to a separate (new) directory
   - the R script (RNA-Seq-final.Rmd)
   - the output from the previous script (GE.Rdata)
   - file with GO pathways for Biological Processes (c5.bp.v7.2.symbols.gmt; the name can be different because it can be updated). This is an universal file downloaded from the GSEA website (https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp). Go to the chapter **C5** standing for GO, then to **BP** subchapter standing for Biological Processes and download the GMT file containing Gene Symbols. 
2. Adjust the R script as needed. Important points are given below. For factor names and possible values, please refer to the report from the previous/initial R script.
   - Load the input file: _GeneExpression <- readRDS("GE.Rdata")_
   - Load the GMT file: _GMTFileName <- "c5.bp.v7.2.symbols.gmt"_ (use correct name if the file was updated)
   - Adjust the list of factors: _FactorsList <- c( "treatment", "sex", "sow", "ab_status")_. Here, provide the list of all factors which you want to include into following analysis or just all factors. The title of factors is coming from the previous R file report.
   - Setting up the GO analysis. Alternatively you can do GO analysis online using the Enrichr (https://maayanlab.cloud/Enrichr/)
     - Select the main factor for GO analysis: e.g. for “treatment” **CurrentFactor <- 1**. The number 1 shows that the main factor for the GO analysis will be #1 in _FactrorsList_ above – “treatment”.
     - Select the FDR limit for GO analysis: _FDRLimit <- 0.1_. This will show genes with FDR of 0.1 or below as significant.
     - Select the favourite GO for barcode and heatmap charts: e.g. _GeneSetOfInterest <- "GO_CEREBELLAR_CORTEX_FORMATION"_
   - You need to select samples you need for the current analysis. As an example the line _GeneExpression <- GeneExpression[, GeneExpression$samples$dev_stage == "juvinile"]_ will preselect only piglets and not tissues. Because only piglets had “juvenile” stage of development. It is mentioned here as “dev_stage” factor. In the next step normalization of selected samples will happen. Consider removing the biased samples (correction factor is more than 5- to 10-fold).
   - First you get unsupervised PCA charts with different colors for different factors. Final PCA will be printed for the same factors as here. You can adjust factors by modifying the FactorsList parameter. Look to the created PCA charts to choose wisely which factor to keep and which ones to delete. To keep or delete factors, make the appropriate model of the study in the next step. Besides this you can consider to exclude some animals (e.g. males or females). This could be done by adding more lines to the previous step. E.g. to exclude females add one of the following lines:
     - GeneExpression <- GeneExpression[, GeneExpression$samples$sex == "male"])
     - GeneExpression <- GeneExpression[, GeneExpression$samples$sex != "female"])
   - Model design. Please refer to the table below for making the model. Please pay attention that for non-binary factors (e.g. for models #4 and #5, Ab status that may have 3 possible outcomes). Factors taken (+) and excluded (-) from the analysis

Model #1. +Infection status (one binary (Yes/No) factor (Treatment) is taken into account. This is the most simple option)
```
ModelDesign <- model.matrix(~treatment, data = GeneExpression$samples)
ExpressionValues <- voom(GeneExpression, ModelDesign, plot = T)
ModelFit <- lmFit(ExpressionValues, ModelDesign) %>% eBayes
Result <- topTable(ModelFit, number = Inf, sort.by = "logFC")
```
Model #2. +Infection status -Gender (two binary factors are taken into account: Treatment and Sex. Model is solved for two factors. But results are given only for one factor. To be more precise for one outcome (_Treatment = infected_). Thus, we exclude effect of the second factor (sex) from the answer.)
```
ModelDesign <- model.matrix(~treatment + sex, data = GeneExpression$samples)
ExpressionValues <- voom(GeneExpression, ModelDesign, plot = T)
ModelFit <- lmFit(ExpressionValues, ModelDesign) %>% eBayes
Result <- topTable(ModelFit, coef = "treatmentinfected", number = Inf, sort.by = "logFC")
```
Model #3. +Infection status with removal of the batch effect for Gender (two binary factors are taken into account: Treatment and Sex. Model is solved. Batch effect of the second factor (Sex) is removed. Results are given only for one factor. To be more precise for one outcome (_Treatment = infected_). Thus, we exclude effect of the second factor (sex). This model is similar to the one above. However, it is doing the same thing with other means.)
```
ModelDesign <- model.matrix(~treatment + sex, data = GeneExpression$samples)
ExpressionValues <- voom(GeneExpression, ModelDesign, plot = T)
ExpressionValues$E <- removeBatchEffect(ExpressionValues, batch = ExpressionValues$targets$sex)
ModelFit <- lmFit(ExpressionValues, ModelDesign) %>% eBayes
Result <- topTable(ModelFit, coef = "treatmentinfected", number = Inf, sort.by = "logFC")
```
Model #4. +Ab status -Gender (This model is similar to #2. But one of the factors (_Ab_status_) is non binary.)
```
ModelDesign <- model.matrix(~ab_status + sex, data = GeneExpression$samples)
ExpressionValues <- voom(GeneExpression, ModelDesign, plot = T)
ModelFit <- lmFit(ExpressionValues, ModelDesign) %>% eBayes
Result <- topTable(ModelFit, coef = c("ab_statusN", "ab_statusS", "ab_statusP"), number = Inf, sort.by = "F")
```
Model #5. +Ab status with removal of the batch effect for Gender (This model is similar to #3. But one of the factors (_Ab_status_) is non binary.)
```
ModelDesign <- model.matrix(~ab_status + sex, data = GeneExpression$samples)
ExpressionValues <- voom(GeneExpression, ModelDesign, plot = T)
ExpressionValues$E <- removeBatchEffect(ExpressionValues, batch = ExpressionValues$targets$sex)
Contrasts <- makeContrasts(PvsNeg = ab_statusP - ab_statusN, SvsNeg = ab_statusS - ab_statusN, levels = colnames(ModelDesign))
ModelFit <- lmFit(ExpressionValues, ModelDesign) %>% eBayes
ModelFit <- contrasts.fit(ModelFit, contrasts = Contrasts) %>% eBayes()
Result <- topTable(ModelFit, number = Inf, sort.by = "F")
```
Notes:
- Models #2-5 are solved for 2 factors, but results are revealed only for one factor. The “coef” parameter is a shortcut for “main factor name + one of the outcomes”. See results of the **print(ModelDesign)** to see possible outcomes and correct names for this variable.
- Try to keep models as simple as possible.
- GO analysis is more sensitive then DEG analysis. You can use more conservative model for GO analysis, even if it gives almost no results for DEGs.

3. Open and run the R script with R-Studio.
   - Use **Code > Run Region > Run All** (Control+Alt+R) if you want to see results in R-Studio.
   - Use **File > Knit Document** if you need a report file to be generated.
4. As a result, several files are created
   - HTML report is created
     - Initial and adjusted PCA charts are included in this report. Raw data for latter one is also in the report.
     - Adjusted gene-expression values (after solving the model) are stored in GE.csv. This file can be used to identify sex in your samples, if necessary. E.g. Teixeira et al. (https://doi.org/10.3390/genes10121010) identified DEGs in porcine fetuses bsed on their sex: DDX3Y, KDM5D, ZFY, EIF2S3Y, EIF1AY. Find them in Excel-file. These genes are expressed well only in males.

![male genes](https://user-images.githubusercontent.com/9166776/155331034-ca80c994-7441-4de3-ba3a-cfec25db7585.png)

Figure 2. Expression of sex-specific marker genes (highlighted with blue bars) in adjusted Gene Expression dataset. In females, their expression is significantly reduced. Expression of several other transcripts (highlighted with red bars) is given as reference to show an example of genes whose expression is independent on sex. 

   - Files with DEGs
     - Data in machine-readable file-format: DEG.Rdata
     - Excel-readable file format of data: DEG.csv. Edit this file with Excel. Export data to other programs for making charts. If you need to make volcano plots, use GraphPad Prism. If you need to make Venn diagrams, use PowerPoint or onlline tools (http://www.interactivenn.net/). If you need to compare lists of DEGs, use https://editor.mergely.com.
   - GO pathways
     - GO.csv with all GOs
     - GO for Cytoscape.txt will contain more detailed information for Cytoscape.

## 7 ENRICHMENT MAP ANALYSIS
Full manual for this application is available online: https://cytoscape.org/cytoscape-tutorials/protocols/enrichmentmap-pipeline/#/.
1. Launch Cytoscape.
2. Create Enrichment Map: **Apps > EnrichmentMaps > “+” button > Select “Generic/gProfiler/Enrichr” for Analysis type > Locate your GO for Cytoscape.txt file for Enrichments option > Put checkmark for Show Advanced Options > Put checkmark for Filter genes b expressions > Type 0.1 for FDR q-value cutoff and 1.0 for p-value cutoff > Check other parameters: Data set edges: automatic, Cutoff: 0.375, Jaccard+Overlap 50:50 > “Build” button.** If you want to have more dense network between clusters, make cut-off lower than 0.375 (e.g. 0.2). If you want less GOs on the chart, make FDR cut-off lower (e.g. 0.01).
3. Adjust formatting: **Select EnrichmentMap Tab on the left side > Select Chart data: NES Columns and Select Publication-Ready**.
4. Delete clusters with <3 elements: **Select Filter Tab on the left side > “+” button > Topology Filter > Select “at least”, “2” neighbours within distance “2”. Menu Select > Nodes > Hide Unselected Nodes. Menu File > New Network > From Selected Nodes, All Edges.**
5. Define clusters of similar pathways representing major biological themes: **Apps > AutoAnnotate > New Annotation Set… > Advanced tab > Select Use clusterMaker App > Select Community cluster (GLay) as Cluster algorithm > Select Layout network to prevent cluster overlap > Type 4 for Max words per label and 8 for Adjacent word bonus > “Create annotations” button > “OK”.**
6. Adjust formatting: **Select AutoAnnotate Display Tab on the right side > Unselect “Scale font by cluster size”, select “Word wrap”.** Select Chart data: **NES Columns** and Select Publication-Ready.
7. Save the image: **Menu File > Export > Network to Image... > Export file format: PDF > Ok.**
 
## 8 SUPPLEMENTARY PARTS (ALTERNATIVE PARTS OF THE PATHWAY)
### 8.1 Trim adaptors and do quality trim with Trimmomatic
This part is optional. You need to run it only in case that you want to start from FASTQ.GZ files containing adaptors.

Install Trimmomatic (http://www.usadellab.org/cms/?page=trimmomatic). Download and unpack it to /home/user/bin/Trimmomatic-0.40/. UGENE contains a build-in version of Trimmomatic and can replace the Linux version.
```
unzip -d ~/bin/ ./Trimmomatic-0.40.zip
```
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
