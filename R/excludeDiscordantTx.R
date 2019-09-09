##### Date: 09 Sep 2019
##### Contact: Trung Nghia Vu (nghiavtr@gmail.com)
##### This script is to remove the discordant transcripts in cdna fasta but not existing in gtf file name from ensembl annotation
##### Command:  Rscript excludeDiscordantTx.R cdna=transcript.fa sqlite=transcript.sqlite out=transcript.clean.fa

args = commandArgs(trailingOnly=TRUE)
for (i in 1:length(args)){
	res=unlist(strsplit(args[i],"="))
	if (res[1]=="cdna") cdnaFn=res[2]
	if (res[1]=="sqlite") sqliteFn=res[2]
	if (res[1]=="out") outFn=res[2]
}


suppressMessages(library("GenomicFeatures"))
suppressMessages(library("Biostrings"))
anntxdb <- loadDb(sqliteFn)
#get tx in gtf
alltx=mcols(transcripts(anntxdb))$tx_name
#read fasta file
fasta_tx = readDNAStringSet(cdnaFn)
#get tx in cdna
cdna_tx=sapply(names(fasta_tx),function(x) unlist(strsplit(x," "))[1])
#get the concordant tx
pick=cdna_tx %in% alltx
exportFasta=fasta_tx[pick]
#export to file
writeXStringSet(exportFasta, outFn)
cat("Done!")
