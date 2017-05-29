#####
# create sqlite from gtf using GenomicFeatures, for example Homo_sapiens.GRCh37.75.sqlite
# Rscript createSqlite.R Homo_sapiens.GRCh37.75.gtf Homo_sapiens.GRCh37.75.sqlite
############################################################
args = commandArgs(trailingOnly=TRUE)
gtfFile=args[1];
gtfSqliteFn=args[2];
library(GenomicFeatures);
gtfTxdb <- makeTxDbFromGFF(file=gtfFile,
                 format="gtf",
                 dataSource=paste("Link to the source",sep=""),
                 organism="Homo sapiens")
saveDb(gtfTxdb,file=gtfSqliteFn)