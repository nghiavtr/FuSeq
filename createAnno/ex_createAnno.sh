##### Create the annotation RData file for Plant: An example for Arabidopsis_thaliana.TAIR10
### Requirements: 
#R packages: org.Hs.eg.db, biomaRt, polyester, Biostrings, foreach and doParallel
#FuSeq: version 1.1.0 or later

### NOTE: Do not forget to replace '/path/to/FuSeq_home' by your correct path
echo "NOTE: Do not forget to replace '/path/to/FuSeq_home' by your correct path"

export LD_LIBRARY_PATH=/path/to/FuSeq_home/lib:$LD_LIBRARY_PATH
export PATH=/path/to/FuSeq_home/bin:$PATH

### 1. Download the cdna fasta file and gtf file of the species
## for example:
wget ftp://ftp.ensembl.org/pub/release-XX/fasta/SPECIES/cdna/SPECIES.version.cdna.all.fa.gz
gunzip SPECIES.version.cdna.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-XX/gtf/SPECIES/cdna/SPECIES.version.gtf.gz
gunzip SPECIES.version.gtf.gz
# normally a fasta file and a gtf file will be created: SPECIES.version.cdna.all.fa and SPECIES.version.gtf

### 2. clean version in transcript names and gene names in Ensembl annotation
# This step need to be done if there are a inconsistency between transcript/gene names in cdna fasta file and gtf file.
Rscript /path/to/FuSeq_home/R/excludeTxVersion.R SPECIES.version.cdna.all.fa SPECIES.version.cdna.all.cleanversion.fa

### 3. build sqlite
Rscript  /path/to/FuSeq_home/R/createSqlite.R SPECIES.version.gtf SPECIES.version.sqlite 

### 4. Remove the discordant transcripts existing in the fasta cdna file but not in the gft file.
Rscript /path/to/FuSeq_home/R/excludeDiscordantTx.R  cdna=SPECIES.version.cdna.all.cleanversion.fa sqlite=SPECIES.version.sqlite  out=SPECIES.version.cdna.all.clean.fa

### 5. run R-script to get a annotation RData file
# PLEASE READ THE NOTES IN ex_createAnno.R FILE BEFORE USE
Rscript ex_createAnno.R SPECIES.version.cdna.all.clean.fa SPECIES.version.c.txAnno.RData

echo "Done!"
