##### Create the annotation RData file for Plant: An example for Arabidopsis_thaliana.TAIR10
### Requirements: 
#R packages: org.Hs.eg.db, biomaRt, polyester, Biostrings, foreach and doParallel
#FuSeq: version 1.1.0 or later

### NOTE: Do not forget to replace '/path/to/FuSeq_home' by your correct path
echo "NOTE: Do not forget to replace '/path/to/FuSeq_home' by your correct path"

export LD_LIBRARY_PATH=/path/to/FuSeq_home/lib:$LD_LIBRARY_PATH
export PATH=/path/to/FuSeq_home/bin:$PATH

# 1. Download cdna fasta file
wget ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz

# 2. clean version in transcript names and gene names
Rscript /path/to/FuSeq_home/R/excludeTxVersion.R Homo_sapiens.GRCh38.cdna.all.fa Homo_sapiens.GRCh38.cdna.all.clean.fa

# 3. run R-script to get a annotation RData file
Rscript ex_createAnnoRData_Homo_sapiens.GRCh38.94.R Homo_sapiens.GRCh38.cdna.all.clean.fa Homo_sapiens.GRCh38.94.txAnno.RData

echo "Done!"

