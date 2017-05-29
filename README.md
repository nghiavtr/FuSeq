#####################################
# Documents for FuSeq version 0.1.0
#####################################

## 1. Introduction
FuSeq is a novel method to discover fusion genes from paired-end RNA sequencing data. We acknowledge for materials from Sailfish, Rapmap and other tools used in this software.

Software requirements for FuSeq:
- A C++-11 compliant compiler version of GCC (g++ >= 4.8.2)
- R packages version 3.3.0 or latter with following installed packages: GenomicFeatures

Annotation reference: FuSeq requires a fasta file of transcript sequences and a gtf file of transcript annotation. FuSeq supports the ensembl annotation version GRCh37.75 in the current version:
- Download the sequences of transcripts: http://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.75.cdna.all.fa.gz
- Download the annotation of transcripts: http://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz

The latest version and information of FuSeq are updated at https://github.com/nghiavtr/FuSeq

## 2. Download and installation
If you use the binary verion of FuSeq: 
- Download the lastest binary version from FuSeq website: [FuSeq_v0.1.0_linux_x86-64](https://github.com/nghiavtr/FuSeq/releases/download/v0.1.0/FuSeq_v0.1.0_linux_x86-64.tar.gz)
- Uncompress to folder
```sh
tar -xzvf FuSeq_v0.1.0_linux_x86-64.tar.gz
```
- Move to the *FuSeq_home* directory and do configuration for FuSeq
```sh
cd FuSeq_v0.1.0_linux_x86-64
bash configure.sh
```
- Add paths of lib folder and bin folder to LD_LIBRARY_PATH and PATH
```sh
export LD_LIBRARY_PATH=/path/to/FuSeq_v0.1.0_linux_x86-64/linux/lib:$LD_LIBRARY_PATH
export PATH=/path/to/FuSeq_v0.1.0_linux_x86-64/linux/bin:$PATH
```
If you want to build FuSeq from sources:
- Download FuSeq from [FuSeq website](https://github.com/nghiavtr/FuSeq) and move to *FuSeq_home* directory
```sh
wget https://github.com/nghiavtr/FuSeq/archive/v0.1.0.tar.gz
tar -xzvf v0.1.0.tar.gz
cd FuSeq-0.1.0
```
- FuSeq requires information of flags from Sailfish including DFETCH_BOOST, DBOOST_ROOT, DTBB_INSTALL_DIR and DCMAKE_INSTALL_PREFIX. Please refer to the Sailfish website for more details of these flags.
- Do installation by the following command:
```sh
DBOOST_ROOT=/path/to/boostDir/ DTBB_INSTALL_DIR=/path/to/tbbDir/ DCMAKE_INSTALL_PREFIX=/path/to/expectedBuildDir bash install.sh
```
-After the installation is finished, remember to add the paths of lib folder and bin folder to LD_LIBRARY_PATH and PATH
```sh
export LD_LIBRARY_PATH=/path/to/expectedBuildDir/lib:$LD_LIBRARY_PATH
export PATH=/path/to/expectedBuildDir/bin:$PATH
```

### Do not forget to replace "/path/to/" by your local path.

## 3. Index
The command for indexing transcript fasta is similar to the indexing of sailfish and rapmap. For example:
```sh
TxIndexer -t /path/to/transcripts.fa -o TxIndexer_idx
```
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
- gene5,chrom5,strand5,gene5.start,gene5.end: information of gene at 5 prime including gene id. chromosome, strand, starting and ending of the fusion junction
- gene3,chrom3,strand3,gene3.start,gene3.end: similar information but for gene at 3 prime
- fusionName: names of fusion genes
- supportRead: the number of reads supporting the fusion genes
- score: score of each fusion genes from FuSeq
- info: additional information of the fusion gene, for example constituent genes involve in mitochondrial translation, cytosolic ribosomal subunit and ribonucleoprotein.

##### - If keepRData=TRUE is set in the *params_file*, all data generated during the FuSeq processes from will be saved in a *.RData file
##### - If exportFasta=TRUE is set in the *params_file*, FuSeq will export fusion reads supporting fusion genes in two fasta files named *_fusionReads_1.fa and *_fusionReads_2.fa for read1 and read2, respectively.

## 6. Parameter setting
There are several parameters input for the pipeline. The default parameter setting is available at FuSeq_home/R/params.txt. Details of these parameters are refered to main publication, herein they are shortly described as the below:
- readStrands (UN): indicates the directions of read1 and read2 depending on library type. There are a total of 5 types FF, FR, RF, RR, UN. "F" and "R" indicate forward and reverse complement directions, respectively. For example, if the libraries are prepared by TruSeq Stranded mRNA LT Sample Prep Kit where read1 is reverse complement and read2 is forward, a setting of readStrands=RF is suitable.The setting of readStrands=UN is used for unstranded sequencing data.
- chromRef: limits the chromosomes contributing to fusion genes, we usually consider chromosomes 1:22,X and Y (chromRef=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y)
- onlyProteinCodingGenes (TRUE): if users want to keep only fusion genes from protein-coding genes.
- maxSharedCount (5e-2): the maximum proportion of sharing counts of a fusion gene.
- minGeneDist (1e5): the minimum distance between two constituent genes.
- minJunctionDist (1e5): the minimum junction distance of fusion gene.
- maxInvertedFusionCount (0.01): the proportion of expression of inverted fusion.
- maxMRfusionFc (2), maxMRfusionNum (2) and sgtMRcount (10): the parameter settings for fitlering multiple-fusion genes in the mapped read (MR) pipeline.
- minMR (2), minNonDupMR (2): the minimum supporting count (non-duplicated supporting count) of a fusion gene in the mapped read pipeline.
- minSR (1): the minimum supporting count of a fusion gene in the split read (SR) pipeline.
- minScore (3): the minimum score of a fusion gene.
- keepRData (FALSE): if users want to save all data generated during fusion gene discovery using Rscripts. If it is set keepRData=TRUE, a single *.RData file will be created to store the data.
- exportFasta (FALSE): Set exportFasta=TRUE if users want to export sequences of fusion reads

## 7. License
FuSeq uses GNU General Public License GPL-3

## 8. References
(update later)
