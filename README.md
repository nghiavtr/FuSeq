#####################################
# Documents for FuSeq version 1.1.3
#####################################

## Update news
#### 19 Aug 2020: version 1.1.3
1) speed up the processSplitRead() in the paralogs checking
2) export fragmentDist.txt for the data with a low mapped read rate
3) fix the error when installing from sources caused by the dependency fmtlib mentioned in https://github.com/nghiavtr/FuSeq/issues/1
4) fix the error when running the binary version in Ubuntu

#### 09 Sep 2019: version 1.1.2
1) Fix small bugs in processSplitRead.R and postProcessSplitRead.R
2) Add excludeDiscordantTx.R (details in Section 7) to remove the discordant transcripts which are existing in fasta cdna file but not in gtf file. These discordant transcripts can crash down the split-read pipeline. This step is recommmended for any untested or new annotations such as Homo_sapiens.GRCh38.

#### 13 June 2019: version 1.1.1
1) Speed up functions processFEQ() and detectJunctionBreaks()
2) Add a parameter (exonBoundary) to control junction breaks inside exons
3) Improve clean codes in FuSeq.cpp

#### 14 Nov 2018: version 1.1.0
1) Add the description of how to generate the annotation files for new species
2) Minor changes in output data of FuSeq

#### 01 Oct 2018: version 1.0.0
1) Correct the headers of split reads when exporting their sequences to file

#### 23 Jul 2018: version 0.1.1
1) Error when running with Ubuntu 16 from librt.so.1: fixed by simply removing the librt.so.1 from the library folder of FuSeq

2) Allow inverted fusion genes by default: there are existing both fusion genes A-B and B-A
maxInvertedFusionCount=0

3) Optimize memory and space
processFEQ.R: remove myFusion as a redundant data
FuSeq_functions.R: fix some bugs of memory usage in function computeGeneDistance()

4) Export more information in the output fusions.FuSeq

#### 21 May 2017: version 0.1.0
- First submission

## 1. Introduction
FuSeq is a novel method to discover fusion genes from paired-end RNA sequencing data. FuSeq discovers fusion genes based on quasi-mapping to quickly map the reads, extract initial candidates from split reads and fusion equivalence classes of mapped reads, and finally apply multiple filters and statistical tests to get the final candidates. Full details of the method is described in the method publication (1).

### Software requirements:
FuSeq is implemented in R and C++. We acknowledge for materials from Sailfish, Rapmap and other tools used in this software.
- A C++-11 compliant compiler version of GCC (g++ >= 4.8.2)
- R packages version 3.3.0 or latter with following installed packages: GenomicFeatures

### Annotation reference
FuSeq requires

1) a fasta file of transcript sequences and a gtf file of transcript annotation: can be downloaded from public repositories such as Ensembl (ensembl.org)
2) a RData file of supporting annotation: A description of how to create the RData file for new annotation versions or species is available in Section 7. 

Current FuSeq version was tested on the human transcriptome with ensembl annotation version GRCh37.75. Following files are required:
- Sequences of transcripts (ensembl website) [GRCh37.75 cdna fasta](http://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.75.cdna.all.fa.gz)
- Gtf annotation of transcripts (ensembl website) [GRCh37.75 gtf annotation](http://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz)
- RData annotation file (FuSeq website): [GRCh37.75 RData annotation](https://github.com/nghiavtr/FuSeq/releases/download/v0.1.0/Homo_sapiens.GRCh37.75.txAnno.RData)

We also prepare annotation RData files for several annotations from Ensembl (ensembl.org) and example codes to build the annotation RData file in folder [*createAnno*](https://github.com/nghiavtr/FuSeq/tree/master/createAnno). **Please see Section 7 for further information.**
- [Homo sapiens version GRCh38.94](https://github.com/nghiavtr/FuSeq/releases/download/v1.1.0/Homo_sapiens.GRCh38.94.txAnno.RData): *However, before running FuSeq for this annotation, users should clean the version information from the transcript names and gene names as described in Section 7*.
- [Arabidopsis thaliana version TAIR10.41](https://github.com/nghiavtr/FuSeq/releases/download/v1.1.0/Arabidopsis_thaliana.TAIR10.41.txAnno.RData)

### Versions
The latest version and information of FuSeq are updated at https://github.com/nghiavtr/FuSeq

The older versions can be found here:
- Version 1.1.2: https://github.com/nghiavtr/FuSeq/releases/tag/v1.1.2
- Version 1.1.1: https://github.com/nghiavtr/FuSeq/releases/tag/v1.1.1
- Version 1.1.0: https://github.com/nghiavtr/FuSeq/releases/tag/v1.1.0
- Version 1.0.0: https://github.com/nghiavtr/FuSeq/releases/tag/v1.0.0
- Version 0.1.1: https://github.com/nghiavtr/FuSeq/releases/tag/v0.1.1
- Version 0.1.0: https://github.com/nghiavtr/FuSeq/releases/tag/v0.1.0

## 2. Download and installation
If you use the binary verion of FuSeq: 
- Download the lastest binary version from FuSeq website: [FuSeq_v1.1.3_linux_x86-64](https://github.com/nghiavtr/FuSeq/releases/download/v1.1.3/FuSeq_v1.1.3_linux_x86-64.tar.gz)
```sh
wget https://github.com/nghiavtr/FuSeq/releases/download/v1.1.3/FuSeq_v1.1.3_linux_x86-64.tar.gz -O FuSeq_v1.1.3_linux_x86-64.tar.gz
```
- Uncompress to folder
```sh
tar -xzvf FuSeq_v1.1.3_linux_x86-64.tar.gz
```
- Move to the *FuSeq_home* directory and do configuration for FuSeq
```sh
cd FuSeq_v1.1.3_linux_x86-64
bash configure.sh
```
- Add paths of lib folder and bin folder to LD_LIBRARY_PATH and PATH
```sh
export LD_LIBRARY_PATH=/path/to/FuSeq_v1.1.3_linux_x86-64/linux/lib:$LD_LIBRARY_PATH
export PATH=/path/to/FuSeq_v1.1.3_linux_x86-64/linux/bin:$PATH
```
If you want to build FuSeq from sources:
- Download FuSeq from [FuSeq website](https://github.com/nghiavtr/FuSeq) and move to *FuSeq_home* directory
```sh
wget https://github.com/nghiavtr/FuSeq/archive/v1.1.3.tar.gz
tar -xzvf v1.1.3.tar.gz
cd FuSeq-1.1.3
```
- FuSeq requires information of flags from Sailfish including DFETCH_BOOST, DBOOST_ROOT, DTBB_INSTALL_DIR and DCMAKE_INSTALL_PREFIX. Please refer to the Sailfish website for more details of these flags.
- Do installation by the following command:
```sh
DBOOST_ROOT=/path/to/boostDir/ DTBB_INSTALL_DIR=/path/to/tbbDir/ DCMAKE_INSTALL_PREFIX=/path/to/FuSeq_home bash install.sh
```
-After the installation is finished, remember to add the paths of lib folder and bin folder to LD_LIBRARY_PATH and PATH
```sh
export LD_LIBRARY_PATH=/path/to/FuSeq_home/lib:$LD_LIBRARY_PATH
export PATH=/path/to/FuSeq_home/bin:$PATH
```

### Do not forget to replace "/path/to/" by your local path.

## 3. Index
The command for indexing transcript fasta is similar to the indexing of sailfish and rapmap. For example:
```sh
TxIndexer -t /path/to/transcripts.fa -o TxIndexer_idx
```

**If the transcript names in the transcript fasta file contains version data but the gtf file does not, users should remove the version information from the transcript names as described in Section 7**

#### K-mer length for a short-read RNA-seq
In order to be able to capture split reads, the k-mer length must be less than a half of read length. Thus, the default k-mer length (31) of rapmap/sailfish is not proper for a short-read RNA-seq dataset, e.g 50 bp long. For the short-read RNA-seq dataset, we use k-mer length of 21 instead as the follow:
```sh
TxIndexer -t /path/to/transcripts.fa -k 21 -o TxIndexer_idx
```


## 4. Extract fusion equivalence classes and split reads
This step requires a annotation file (gtf_file) in the gtf format to get the map between gene and transcript.
- The command to generate  fusion equivalence classes is similar to "sailfish quant". For example, we want to run FuSeq with 8 cpus:
```sh
FuSeq -i /path/to/TxIndexer_idx -l IU -1 read1.fasta -2 read2.fasta -p 8 -g /path/to/gtf_file -o /path/to/feqDir
```
- If the data is compressed in gz format. We can combine with gunzip for depressing on-fly:
```sh
FuSeq -i /path/to/TxIndexer_idx -l IU -1 <(gunzip -c read1.gz) -2 <(gunzip -c read2.gz) -p 8 -g /path/to/gtf_file -o /path/to/feqDir
```
The library type (-l option) in these command is inherited from the [saifish/salmon tools](https://salmon.readthedocs.io/en/latest/library_type.html). The default value "IU" indicates for the unstranded paired-end library which is suitable to most cases.

## 5. Discover fusion genes
FuSeq uses R scripts in *FuSeq_home/R* directory for this step. The transcript fasta (/path/to/transcripts.fa in the index step) and two additional files are required:
- A sqlite file containing the annotation file created by using function makeTxDbFromGFF() of GenomicFeatures. We prepare simple R codes (createSqlite.R) to generate the sqlite file, for example: 
```sh
Rscript FuSeq_home/R/createSqlite.R Homo_sapiens.GRCh37.75.gtf Homo_sapiens.GRCh37.75.sqlite 
```
- A supporting annotation file containing information of paralogs, gene types, etc. For human Ensemble annotation version  GRCh37.75, the file was prepared and available for download ([Homo_sapiens.GRCh37.75.txAnno.RData](https://github.com/nghiavtr/FuSeq/releases/download/v0.1.0/Homo_sapiens.GRCh37.75.txAnno.RData)). The similar files for other annotation versions are preparing and will be available soon.
```sh
 wget https://github.com/nghiavtr/FuSeq/releases/download/v0.1.0/Homo_sapiens.GRCh37.75.txAnno.RData -O Homo_sapiens.GRCh37.75.txAnno.RData
```

Now, we can use another R script to discover fusion genes.
```sh
Rscript FuSeq_home/R/FuSeq.R in=/path/to/feqDir txfasta=/path/to/transcripts.fa sqlite=/path/to/sqlite_file txanno=/path/to/txanno_file out=FuSeq_outDir params=/path/to/params_file
```
Herein, the settings for "params" and "out" are optional. If the params_file is not supplied, the default parameter setting will be used. Details about parameter setting are presented in next section. If FuSeq_outDir is not input, FuSeq will save the output into the same directory of the fusion equivalence classes (/path/to/feqDir).
Since this step is usually fast, parallel is not implemented.

### Output of FuSeq
Fusion genes discovered by FuSeq are stored in a file named *fusions.fuseq* containing information of the fusion genes:
- gene5, chrom5, strand5, brpos5, and cds.brpos5.start: information of gene at 5 prime including gene id, chromosome, strand, the exon-boundary breaking-point of the 5 prime side and the breaking-point of the 5 prime side derived from mapped transcript sequence. The position of cds.brpos5.start is obtained directly from read-mapping information, however it is sometimes get shifted a bit from the true breaking points due to sequence similarities (in SR) or lack of reads nearby (in MR). Usage of cds.brpos5.start is might be interesting for investigate special breaking points, for instance inside exon sequence.
- gene3, chrom3, strand3, brpos3, and cds.brpos3.start: similar information but for the gene at 3 prime
- fusionName: names of fusion genes
- symbol5 and symbol3: the gene symbol of gene5 and gene3, respectively
- SR.passed and MR.passed: the number of reads passed the SR and MR pipeline, respectively
- supportRead: the number of reads supporting the fusion genes
- score: score of each fusion genes from FuSeq
- info: if available, additional information of the fusion gene is reported here. For example whether the constituent genes involve in mitochondrial translation, cytosolic ribosomal subunit and ribonucleoprotein.

A log file (FuSeq_logs.RData) contains information of parameter (*FuSeq.params*), input data path (*inPath*) and the results of final integration step (*FuSeq.integration*). *FuSeq.integration* contains the final results of MR pipeline, SR pipeline and combination.

##### - If exportFasta=TRUE is set in the *params_file*, FuSeq will export fusion reads supporting fusion genes in two fasta files named \*_fusionReads_1.fa and \*_fusionReads_2.fa for read1 and read2, respectively.

##### - If keepRData=TRUE is set in the *params_file*, FuSeq will store results of each steps of the pipeline in \*.RData files.

## 6. Parameter setting
There are several parameters input for the pipeline. The default parameter setting is available at FuSeq_home/R/params.txt. Details of these parameters are refered to main publication, herein they are shortly described as the below:
- readStrands (UN): indicates the directions of read1 and read2 depending on library type. There are a total of 5 types FF, FR, RF, RR, UN. "F" and "R" indicate forward and reverse complement directions, respectively. For example, if the libraries are prepared by TruSeq Stranded mRNA LT Sample Prep Kit where read1 is reverse complement and read2 is forward, a setting of readStrands=RF is suitable.The setting of readStrands=UN is used for unstranded sequencing data.
- chromRef: limits the chromosomes contributing to fusion genes, we usually consider chromosomes 1:22,X and Y (chromRef=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y)
- onlyProteinCodingGenes (TRUE): if users want to keep only fusion genes from protein-coding genes.
- maxSharedCount (5e-2): the maximum proportion of sharing counts of a fusion gene.
- minGeneDist (1e5): the minimum distance between two constituent genes.
- minJunctionDist (1e5): the minimum junction distance of fusion gene.
- maxInvertedFusionCount (0): the proportion of expression of inverted fusion. If this value is zero (maxInvertedFusionCount=0 by default), we accept both fusion gene A-B and its inverted fusion gene B-A. If, for example, we set maxInvertedFusionCount=0.01, we accept the read count of B-A is at most 1% of the read count of A-B.
- maxMRfusionFc (2), maxMRfusionNum (2) and sgtMRcount (10): the parameter settings for fitlering multiple-fusion genes in the mapped read (MR) pipeline.
- minMR (2), minNonDupMR (2): the minimum supporting count (non-duplicated supporting count) of a fusion gene in the mapped read pipeline.
- minSR (1): the minimum supporting count of a fusion gene in the split read (SR) pipeline.
- minScore (3): the minimum score of a fusion gene.
- exonBoundary (TRUE): Set exonBoundary=FALSE if users allow fusions with breaking point inside exon regions.
- keepRData (TRUE): Set keepRData=TRUE (default) if users want to save all data generated during fusion gene discovery using Rscripts, then several *.RData files contain processed data will be exported.
- exportFasta (FALSE): Set exportFasta=TRUE if users want to export sequences of fusion reads


## 7. Create a supporting annotation RData file
In this section, we describe details how to create a annotation RData for your own annotation. We use the annotations from Ensembl (ensembl.org) as example.

#### Requirements: 
- R packages: org.Hs.eg.db, biomaRt, polyester, Biostrings, foreach and doParallel
- FuSeq: version 1.1.0 or later


#### Concordances of transcript/gene names in the fasta file and the gft file of Ensembl annotation
- *Fix the inconsistency between transcript names and gene names in your fasta cdna file and gft file.*

In some recent annotations from Ensembl (e.g., GRCh38 of human), versions are also included in the transcript and gene names of the fasta file (cdna) of transcript sequences. However, it is not consistent with the corresponding gtf file where there are no version information in the names. To solve this issue, we exclude the version information in transcript names and gene names from the cdna fasta file (transcripts.fa) using excludeTxVersion.R:

```sh
Rscript FuSeq_home/R/excludeTxVersion.R /path/to/transcripts.fa /path/to/transcripts.cleanversion.fa
```

- *Remove the discordant transcripts existing in the fasta cdna file but not in the gft file.*

Furthermore, there might have also discordant transcripts that they are existing in the fasta cdna file but not the gtf file. We remove these discordant transcripts from the fasta cdna file using excludeDiscordantTx.R:

```sh
#creat the sqlite file
Rscript FuSeq_home/R/createSqlite.R transcript.gtf transcript.sqlite 
#remove discordant transcripts from the cdna fasta file
Rscript FuSeq_home/R/excludeDiscordantTx.R  cdna=/path/to/transcripts.cleanversion.fa sqlite=transcript.sqlite out=/path/to/transcripts.clean.fa
```

The fasta cdna file after removing the version information and discordant transcripts (transcripts.clean.fa) should be used for both indexing step and/or creating RData annotation file in downstream.

#### Generation of an annotation RData file
We use the transcripts.fa (or the transcripts.clean.fa from the previous step, if the transcript/gene names are not consistent) to generate the annotation RData file. This includes three main steps:

1) extract information of transcripts such as transcript name, gene name, chromosom etc.
- The information is extracted from the head of sequences in the clean cdna fasta file (transcripts.clean.fa).
- The data are kept in two tables of *txAnno* for transcripts and *geneAnno* for non-duplicated genes. Their columns includes (in the same order) transcript name, cdna type,chromosome name, gene biotype, transcript biotype and gene name.

2) using biomart to get gene paralog and other information
The information of gene paralog is compulsory, all the other information is optional. The data are kept in the objects of: 
- geneParalog (*compulsory*): information of gene paralog
- hgncName (optional): matched gene symbols
- ribSubunitDb (optional): cytosolic ribosomal subunit
- mitoTransDb (optional): mitochondrial translation
- ribonuproDb (optional): ribonucleoprotein complex

3) simulate data and find relations between transcripts
For simplicity, in this manual, we use the R-package "polyester" to simulate data of all transcripts, then run GenTC to generate data for discovering sequence similarities between two genes. 

#### Some remarks in Step 2 before use:
- check to use a correct dataset and version from ensembl database
- check to use a correct attribute containing gene paralog when using biomaRt
- skip the codes for hgncName, ribSubunitDb, mitoTransDb and ribonuproDb if their information is not available

#### Example codes for creating an annotation RData
We also prepare example codes in folder *createAnno* to create a new supporting annotation RData file for 

1) the annotation from an animal (Homo sapiens version GRCh38.94) with inconsistent transcript/gene names between fasta file and gtf file
2) the annotation from a plant (Arabidopsis thaliana version TAIR10.41) with consistent transcript/gene names.

Note that "polyester" R-package might require a huge of memory for generating the simulated data for big transcriptome. In our implementation, we assign a memory of 32GB for TAIR10.41 and 64GB for GRCh38.94.

## 8. Pratical examples of using FuSeq

In this section, we introduce how to use FuSeq by a step-by-step tutorial using several public real RNA-seq samples. This aims to test FuSeq by just doing a copy-and-paste of the example commands.
For simplicity, in this practice, the FuSeq software, the annotation, RNA-seq data samples, and the results of fusion detection are put togetther in the same (current) folder.

### 8.1. Download and install
#### Download and configure FuSeq
```sh
wget https://github.com/nghiavtr/FuSeq/releases/download/v1.1.3/FuSeq_v1.1.3_linux_x86-64.tar.gz -O FuSeq_v1.1.3_linux_x86-64.tar.gz
tar -xzvf FuSeq_v1.1.3_linux_x86-64.tar.gz
cd FuSeq_v1.1.3_linux_x86-64
bash configure.sh
cd ..
```
#### Set paths to FuSeq
```sh
export LD_LIBRARY_PATH=$PWD/FuSeq_v1.1.3_linux_x86-64/linux/lib:$LD_LIBRARY_PATH
export PATH=$PWD/FuSeq_v1.1.3_linux_x86-64/linux/bin:$PATH
```
### 8.2. Download and prepare the reference files
#### Download the fasta and gtf of transcripts
```sh
wget http://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.75.cdna.all.fa.gz
gunzip Homo_sapiens.GRCh37.75.cdna.all.fa.gz
wget http://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
gunzip Homo_sapiens.GRCh37.75.gtf.gz
```
#### Create sqlite
```sh
Rscript FuSeq_v1.1.3_linux_x86-64/R/createSqlite.R Homo_sapiens.GRCh37.75.gtf Homo_sapiens.GRCh37.75.sqlite 
```
#### Download the extra transcript information and annotation from FuSeq
```sh
wget https://github.com/nghiavtr/FuSeq/releases/download/v0.1.0/Homo_sapiens.GRCh37.75.txAnno.RData
```
### 8.3. Parameter setting
The default of parameter setting is located at FuSeq_v1.1.3_linux_x86-64/R/params.txt that we will use for the pratical examples.
For more advanced-level users, we suggest running FuSeq with the setting of keepRData=TRUE to keep the processed data of FuSeq, then FuSeq will save all data into file FuSeq_process.RData. This file contains the results of both mapped read pipeline and split read pipeline, and extra relevant information of fusion gene candidates such as supporting exons, read mapping positions, sequence reads, etc.

### 8.4. An example for a short read sample 
We select a breast cancer sample from the KPL-4 cell line. This sample contains short reads (50bp), a small library size (6.8 M read pairs) that is suitable for testing purpose. For this dataset, to generate split reads, the default k-mer length of 31 is not appropriate. We use k-mer length of 21 instead.
For data downloading and fusion gene detection, it takes around 0.5 hour in total using a single cpu. The time can vary due to download speeds as well.
#### Download dataset
```sh
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR064/SRR064287/SRR064287_1.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR064/SRR064287/SRR064287_2.fastq.gz
```
#### Indexing fasta sequences
Note: because of short read, a small k-mer length k=21 (k=31 by default) is used for this dataset. The index reference can be reused for other short read samples
```sh
TxIndexer -t Homo_sapiens.GRCh37.75.cdna.all.fa -o TxIndexer_idx_k21 -k 21
```
#### Extract fusion equivalence classes and split reads
```sh
FuSeq -i TxIndexer_idx_k21 -l IU -1 <(gunzip -c SRR064287_1.fastq.gz) -2 <(gunzip -c SRR064287_2.fastq.gz) -p 1 -g Homo_sapiens.GRCh37.75.gtf -o SRR064287_feqDir 
```
#### Discover fusion genes
```sh
Rscript FuSeq_v1.1.3_linux_x86-64/R/FuSeq.R in=SRR064287_feqDir txfasta=Homo_sapiens.GRCh37.75.cdna.all.fa sqlite=Homo_sapiens.GRCh37.75.sqlite txanno=Homo_sapiens.GRCh37.75.txAnno.RData out=SRR064287_FuseqOut params=FuSeq_v1.1.3_linux_x86-64/R/params.txt
```
The results is a list of gene-fusion candidates stored in file fusions.FuSeq in the output folder (SRR064287_FuseqOut).

### 8.5. An example for a long read sample 
Now we apply FuSeq for a long read sample from a glioma dataset.
The selected sample SRR934746 contains 24.4 M read pairs with 100bp read long.
For this dataset, the default k-mer length of 31 is used.
It takes around 01-02 hours for both downloading data and detecting fusion gene candidates using a single cpu.
The time also varies depending on download speeds.
#### Download dataset
```sh
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR934/SRR934746/SRR934746_1.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/SRR934/SRR934746/SRR934746_2.fastq.gz
```
#### Indexing fasta sequences
The default k-mer length (k=31) is used for reference indexing.
The index reference can be reused for other long read samples
```sh
TxIndexer -t Homo_sapiens.GRCh37.75.cdna.all.fa -o TxIndexer_idx_k31 -k 31
```
#### Extract fusion equivalence classes and split reads
```sh
FuSeq -i TxIndexer_idx_k31 -l IU -1 <(gunzip -c SRR934746_1.fastq.gz) -2 <(gunzip -c SRR934746_2.fastq.gz) -p 1 -g Homo_sapiens.GRCh37.75.gtf -o SRR934746_feqDir 
```
#### Discover fusion genes
```sh
Rscript FuSeq_v1.1.3_linux_x86-64/R/FuSeq.R in=SRR934746_feqDir txfasta=Homo_sapiens.GRCh37.75.cdna.all.fa sqlite=Homo_sapiens.GRCh37.75.sqlite txanno=Homo_sapiens.GRCh37.75.txAnno.RData out=SRR934746_FuseqOut params=FuSeq_v1.1.3_linux_x86-64/R/params.txt
```
The results is a list of fusion gene candidates stored in file fusions.FuSeq in the output folder (SRR934746_FuseqOut).
## 9. License
FuSeq uses GNU General Public License GPL-3

## 10. References
1. Vu, Trung Nghia, Wenjiang Deng, Quang Thinh Trac, Stefano Calza, Woochang Hwang, and Yudi Pawitan. 2018. “A Fast Detection of Fusion Genes from Paired-End RNA-Seq Data.” BMC Genomics 19 (1): 786. https://doi.org/10.1186/s12864-018-5156-1.


