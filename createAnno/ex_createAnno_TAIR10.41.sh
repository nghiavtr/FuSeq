##### Create the annotation RData file for Plant: An example for Arabidopsis_thaliana.TAIR10
### Requirements: 
#R packages: org.Hs.eg.db, biomaRt, polyester, Biostrings, foreach and doParallel
#FuSeq: version 1.1.0 or later

### NOTE: Do not forget to replace '/path/to/FuSeq_home' by your correct path
echo "NOTE: Do not forget to replace '/path/to/FuSeq_home' by your correct path"

export LD_LIBRARY_PATH=/path/to/FuSeq_home/lib:$LD_LIBRARY_PATH
export PATH=/path/to/FuSeq_home/bin:$PATH

# 1. Download cdna fasta file
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-41/fasta/arabidopsis_thaliana/cdna/Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz
gunzip Arabidopsis_thaliana.TAIR10.cdna.all.fa.gz

## 2. clean version in transcript names and gene names in Ensembl annotation
# This step need to be done if there are a inconsistency between transcript/gene names in cdna fasta file and gtf file. Homo_sapiens.GRCh38.94 is a typical example.
# However, in Arabidopsis_thaliana.TAIR10, the names in cdna fasta file and gtf file are consistent, so this step can be skipped
#Rscript /path/to/FuSeq_home/R/excludeTxVersion.R Arabidopsis_thaliana.TAIR10.cdna.all.fa Arabidopsis_thaliana.TAIR10.cdna.all.clean.fa

# 3. run R-script to get a annotation RData file
Rscript ex_createAnno_TAIR10.41.R Arabidopsis_thaliana.TAIR10.cdna.all.fa Arabidopsis_thaliana.TAIR10.41.txAnno.RData

echo "Done!"
