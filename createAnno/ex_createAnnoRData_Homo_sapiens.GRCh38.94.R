##### Date: 08 Nov 2018
##### Contact: Trung Nghia Vu (nghiavtr@gmail.com)
##### Create the annotation RData file for Animal: An example for Human GRCh38.94
### Requirements: 
#R packages: org.Hs.eg.db, biomaRt, polyester, Biostrings, foreach and doParallel
#FuSeq: version 1.1.0 or later

##### Some remarks in Step 2 before use:
# 1) check to use the correct dataset and version from ensembl database
# 2) check to use the correct attribute containing gene paralog when using biomaRt
# 3) Step 2.2 can be IGNORED for some species that do not contain information of hgncName, ribSubunitDb, mitoTransDb, ribonuproDb

args = commandArgs(trailingOnly=TRUE)
cdnaFn=args[1]
txAnnoFn=args[2]
#cdnaFn="Homo_sapiens.GRCh38.cdna.all.fa"
#cdnaFn="Homo_sapiens.GRCh38.cdna.all.clean.fa"
#txAnnoFn="Homo_sapiens.GRCh38.94.txAnno.RData"
cat("\n Input fasta: ",cdnaFn)
cat("\n Output: ",txAnnoFn)

##### Step 1) get transcript/gene information from the fasta file
con <- file(cdnaFn, "r", blocking = FALSE)
mydata=readLines(con)
close(con)

txStart=unlist(lapply(mydata, function(x) startsWith(x,">")))
txStartID=which(txStart)
txInfo=mydata[txStartID]
#parsing
chrInfo=unlist(lapply(txInfo, function(x) { 
  y=unlist(strsplit(x," "))
  z=unlist(strsplit(y[3],":"))
  return(z[3])
  }))

cdnaType=unlist(lapply(txInfo, function(x) { 
  y=unlist(strsplit(x," "))
  z=unlist(strsplit(y[grep("cdna",y)],":"))
  return(z[length(z)])
  }))

geneType=unlist(lapply(txInfo, function(x) { 
  y=unlist(strsplit(x," "))
  z=unlist(strsplit(y[grep("gene_biotype",y)],":"))
  return(z[length(z)])
  }))

txType=unlist(lapply(txInfo, function(x) { 
  y=unlist(strsplit(x," "))
  z=unlist(strsplit(y[grep("transcript_biotype",y)],":"))
  return(z[length(z)])
  }))

txName=unlist(lapply(txInfo, function(x) { 
  y=unlist(strsplit(x," "))
  return(substring(y[1],2))
}))

geneName=unlist(lapply(txInfo, function(x) { 
  y=unlist(strsplit(x," "))
  z=unlist(strsplit(y[grep("^gene:",y)],":"))
  return(z[length(z)])
}))

#remove version of gene/transcript
txName=unlist(lapply(txName, function(x) unlist(strsplit(x,"\\."))[1]))
geneName=unlist(lapply(geneName, function(x) unlist(strsplit(x,"\\."))[1]))

#put them together
txAnno=cbind(txName,cdnaType,chrInfo,geneType,txType,geneName)
geneAnno=txAnno
dupID=duplicated(as.character(geneAnno[,6]))
geneAnno=geneAnno[!dupID,]
dim(geneAnno)

##### Step 2) using biomart to get gene paralog and other information
#Note: the annotation version (94) must be correct
### Step 2.1) get paralog
# NOTE: 
# NEED TO CHECK THE ATTRIBUTE in biomart data that containing GENE PARALOG INFORMATION
# USUALLY it is with the name format as "SPECIES_paralog_ensembl_gene", for example hsapiens_paralog_ensembl_gene, athaliana_eg_paralog_ensembl_gene, etc
library(biomaRt)
#remove version in the gene name before querying
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",version = 94) #GRCh=38
geneParalog <- getBM(attributes=c('ensembl_gene_id',"hsapiens_paralog_ensembl_gene"), filters = 'ensembl_gene_id', values = geneAnno[,6], mart = ensembl)
pick=!is.na(geneParalog[,2])
geneParalog=geneParalog[pick,]
pick=geneParalog[,2]!=""
geneParalog=geneParalog[pick,]


### Step 2.2) get extra information
# Optional: can be ignored for some species that do not contain information of hgncName, ribSubunitDb, mitoTransDb, ribonuproDb
hgncName <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol', "chromosome_name"),filters = 'ensembl_gene_id', values = geneAnno[,6], mart = ensembl)
library(org.Hs.eg.db)
xx <- as.list(org.Hs.egGO2ALLEGS)
#cytosolic ribosomal subunit  http://www.ebi.ac.uk/QuickGO/GTerm?id=GO:0022627
largeSubunit<-unlist(xx['GO:0022625'])#cytosolic large ribosomal subunit
smallSubunit<-unlist(xx['GO:0022627'])#cytosolic small ribosomal subunit
myGenes=c(as.character(smallSubunit),as.character(largeSubunit))
res <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol', "chromosome_name","entrezgene"),filters = 'entrezgene', values = myGenes, mart = ensembl)
#ribSubunitDb=res[res$chromosome_name%in%chromRef,]
ribSubunitDb=res
#mitochondrial translation
mitoTrans<-unlist(xx['GO:0032543'])
myGenes=as.character(mitoTrans)
res <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol', "chromosome_name","entrezgene"),filters = 'entrezgene', values = myGenes, mart = ensembl)
#mitoTransDb=res[res$chromosome_name%in%chromRef,]
mitoTransDb=res
#ribonucleoprotein complex
ribonupro<-unlist(xx['GO:1990904'])
myGenes=as.character(ribonupro)
res <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol', "chromosome_name","entrezgene"),filters = 'entrezgene', values = myGenes, mart = ensembl)
#ribonuproDb=res[res$chromosome_name%in%chromRef,]
ribonuproDb=res


##### Step 3: now simulate data and find relations between transcripts
### Step 3.1: data simulation
library(polyester)
library(Biostrings)
simoutdir="simulationData"
fasta = readDNAStringSet(cdnaFn)
unif.countmat=width(fasta)*2
#unif.countmat=ifelse(unif.countmat<1000,1000,unif.countmat)
unif.countmat=as.matrix(unif.countmat)
# simulate reads:
simulate_experiment_countmat(cdnaFn, readmat=unif.countmat, outdir=simoutdir,error_rate=0.0, strand_specific=FALSE) 
# we can be specific to the information of read length, fragment length, see manual of Polyester
#simulate_experiment_countmat(cdnaFn, readmat=unif.countmat, readlen=100, fraglen=250, fragsd=50, outdir='simulationData',error_rate=0.0, strand_specific=FALSE) 

### Step 3.2: Now generate equivalence classes
#register the number of CPU
library(foreach)
library(doParallel)
ncores = detectCores()
CPUSNUM = ncores
registerDoParallel(cores=CPUSNUM)
#do indexing
INDEXDIR="idxDir"
cmd_idx1=paste("TxIndexer -t ",cdnaFn," -o ",INDEXDIR, " -f ",sep="")
system(cmd_idx1)
#run GenTC
read1=paste(simoutdir,"/sample_01_1.fasta",sep="")
read2=paste(simoutdir,"/sample_01_2.fasta",sep="")
outGenTCdir="eqcDir"
cmd_GenTC=paste("GenTC -i ",INDEXDIR," -l IU -1 ",read1," -2 ",read2, " -o ",outGenTCdir," -p ",CPUSNUM,sep="")
system(cmd_GenTC)

#### sequence similarity between two genes
txeqFn=paste(outGenTCdir,"/eqClass.txt",sep="")
txeq=read.csv(txeqFn, header =TRUE, sep="\t")
txeq$geneID=txAnno[match(as.character(txeq$Transcript),txAnno[,1]),6]
length(unique(txeq$eqClass))
geeq=txeq[,c(6:7,3)]
myDup=duplicated(geeq)
geeq=geeq[!myDup,]
dim(geeq)
#index from ge to eq
eqgeMap=tapply(geeq[,2],geeq[,1],c)
eqgeMap.len=lapply(eqgeMap,function(x) length(x[which(!is.na(x))]))
keepID=which(eqgeMap.len>1)
length(keepID)
eqgeMap=eqgeMap[keepID]
eqgeMap.len=eqgeMap.len[keepID]
eqgeMap=lapply(eqgeMap,function(x) x[which(!is.na(x))])
#index from eq to ge
geeqMap=tapply(geeq[,1],geeq[,2],c)
geeqMap=lapply(geeqMap,function(x) as.integer(unlist(x)))
length(geeqMap) 

##### Step 4: save to file
save(geeq,eqgeMap,geeqMap,txAnno,geneAnno,geneParalog,hgncName,ribSubunitDb,mitoTransDb,ribonuproDb,file=txAnnoFn)

