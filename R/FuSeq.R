##############
##### analyze input: Rscript test.R in=inputFeqDir txfasta=txFastaFile sqlite=gtfSqlite txanno=txAnnofile  out=outputDir params=paramsFn
inputFeqDir=txFastaFile=gtfSqlite=txAnnofile=outputDir=paramsFn=NA
#get command information
args = commandArgs(trailingOnly=TRUE)
cat("\nNumber of arguments: ",length(args))
cat("\nList of arguments: ",args)
for (i in 1:length(args)){
	res=unlist(strsplit(args[i],"="))
	if (res[1]=="in") inputFeqDir=res[2]
	if (res[1]=="txfasta") txFastaFile=res[2]
	if (res[1]=="sqlite") gtfSqlite=res[2]
	if (res[1]=="txanno") txAnnofile=res[2]
	if (res[1]=="params") paramsFn=res[2]
	if (res[1]=="out") outputDir=res[2]
}



#check input information
validatedCommand=TRUE
if (is.na(inputFeqDir)){
	cat("\nThere is no input folder. Stop!")
	validatedCommand=FALSE
}
if (is.na(txFastaFile)){
  cat("\nThere is no transcript fasta file. Stop!")
  validatedCommand=FALSE
}
if (is.na(gtfSqlite)){
	cat("\nThere is no sqlite file. Stop!")
	validatedCommand=FALSE
}
if (is.na(txAnnofile)){
	cat("\nThere is no txAnno file. Stop!")
	validatedCommand=FALSE
}

if (is.na(outputDir)){
	cat("\n-----")
	cat("\nThere is no output directory from user, the output will be saved into the input directory of fusion equivalence classes.")
	outputDir=inputFeqDir
}else{
	if (dir.exists(outputDir)){
		cat("\nWarning: The output directory is already existed, old results will be over written")
	}else dir.create(outputDir)
}

if (is.na(paramsFn)){
	cat("\n-----")
	cat("\nThere is no params file. Default settings will be used.")
	FuSeq.params=list()
	FuSeq.params$readStrands="UN"
	FuSeq.params$chromRef=as.character(c(1:22,"X","Y"))
	FuSeq.params$onlyProteinCodingGenes=TRUE
	FuSeq.params$maxSharedCount=5e-2
	FuSeq.params$minGeneDist=1e5
	FuSeq.params$minJunctionDist=1e5
	FuSeq.params$maxInvertedFusionCount=0
	FuSeq.params$maxMRfusionFc=2
	FuSeq.params$maxMRfusionNum=2
	FuSeq.params$sgtMRcount=10
	FuSeq.params$minMR=2
	FuSeq.params$minNonDupMR=2
	FuSeq.params$minSR=1
	FuSeq.params$minScore=3
	FuSeq.params$exonBoundary=TRUE
	FuSeq.params$keepRData=TRUE
	FuSeq.params$exportFasta=FALSE
}else{
	paramIn=read.table(paramsFn, sep="=", header=FALSE)
	FuSeq.params=list()
	FuSeq.params$readStrands=as.character(paramIn[which(paramIn[,1]=="readStrands"),2])
	FuSeq.params$chromRef=trimws(unlist(strsplit(as.character(paramIn[which(paramIn[,1]=="chromRef"),2]),",")))
	FuSeq.params$onlyProteinCodingGenes=as.logical(as.character(paramIn[which(paramIn[,1]=="onlyProteinCodingGenes"),2]))
	FuSeq.params$maxSharedCount=as.double(as.character(paramIn[which(paramIn[,1]=="maxSharedCount"),2]))
	FuSeq.params$minGeneDist=as.double(as.character(paramIn[which(paramIn[,1]=="minGeneDist"),2]))
	FuSeq.params$minJunctionDist=as.double(as.character(paramIn[which(paramIn[,1]=="minJunctionDist"),2]))
	FuSeq.params$maxInvertedFusionCount=as.double(as.character(paramIn[which(paramIn[,1]=="maxInvertedFusionCount"),2]))
	FuSeq.params$maxMRfusionFc=as.double(as.character(paramIn[which(paramIn[,1]=="maxMRfusionFc"),2]))
	FuSeq.params$maxMRfusionNum=as.double(as.character(paramIn[which(paramIn[,1]=="maxMRfusionNum"),2]))
	FuSeq.params$sgtMRcount=as.double(as.character(paramIn[which(paramIn[,1]=="sgtMRcount"),2]))
	FuSeq.params$minMR=as.double(as.character(paramIn[which(paramIn[,1]=="minMR"),2]))
	FuSeq.params$minNonDupMR=as.double(as.character(paramIn[which(paramIn[,1]=="minNonDupMR"),2]))
	FuSeq.params$minSR=as.double(as.character(paramIn[which(paramIn[,1]=="minSR"),2]))
	FuSeq.params$minScore=as.double(as.character(paramIn[which(paramIn[,1]=="minScore"),2]))
	FuSeq.params$exonBoundary=as.logical(as.character(paramIn[which(paramIn[,1]=="exonBoundary"),2]))
	FuSeq.params$keepRData=as.logical(as.character(paramIn[which(paramIn[,1]=="keepRData"),2]))
	FuSeq.params$exportFasta=as.logical(as.character(paramIn[which(paramIn[,1]=="exportFasta"),2]))
}



if (validatedCommand){
	cat("\n-----")
	cat("\nParameter settings:")	
	cat("\n readStrands=",FuSeq.params$readStrands)
	cat("\n chromRef="); cat(FuSeq.params$chromRef,sep = ",")
	cat("\n maxSharedCount=",FuSeq.params$maxSharedCount)
	cat("\n onlyProteinCodingGenes=",FuSeq.params$onlyProteinCodingGenes)
	cat("\n minGeneDist=",FuSeq.params$minGeneDist)
	cat("\n minJunctionDist=",FuSeq.params$minJunctionDist)
	cat("\n maxInvertedFusionCount=",FuSeq.params$maxInvertedFusionCount)
	cat("\n maxMRfusionFc=",FuSeq.params$maxMRfusionFc)
	cat("\n maxMRfusionNum=",FuSeq.params$maxMRfusionNum)
	cat("\n sgtMRcount=",FuSeq.params$sgtMRcount)
	cat("\n minMR=",FuSeq.params$minMR)
	cat("\n minNonDupMR=",FuSeq.params$minNonDupMR)
	cat("\n minSR=",FuSeq.params$minSR)
	cat("\n minScore=",FuSeq.params$minScore)
	cat("\n exonBoundary=",FuSeq.params$exonBoundary)
	cat("\n keepRData=",FuSeq.params$keepRData)
	cat("\n exportFasta=",FuSeq.params$exportFasta)	
}

fragDistFn=paste(inputFeqDir,"/fragmentDist.txt",sep="")
if (!file.exists(fragDistFn)) {
	cat("There is no fragmentDist.txt file, your input sequencing data are probably too small. Stop!")
	validatedCommand=FALSE
}

cat("\n-------------------------\n")

if (validatedCommand){
	#load gtf annotation information
	suppressMessages(library("GenomicFeatures"))
	anntxdb <- loadDb(gtfSqlite)
	load(txAnnofile)
	#load R functions
	source("/path/to/FuSeq_functions.R")
	source("/path/to/processFEQ.R")
	source("/path/to/detectJunctionBreaks.R")
	source("/path/to/doBiologicalFilter.R")
	source("/path/to/processMappedRead.R")
	source("/path/to/processSplitRead.R")
	source("/path/to/postProcessMappedRead.R")
	source("/path/to/postProcessSplitRead.R")
	source("/path/to/integrateFusion.R")
	
	FuSeq.params$outputDir=outputDir
	inPath=inputFeqDir
	myFusionOut=NULL;
	
	FuSeq.MR=processMappedRead(inPath,geneAnno=geneAnno,  anntxdb=anntxdb, geeqMap=geeqMap,FuSeq.params=FuSeq.params)
	FuSeq.SR=processSplitRead(inPath,geneAnno=geneAnno, anntxdb=anntxdb, FuSeq.params=FuSeq.params, txFastaFile=txFastaFile)
	
	FuSeq.MR.postPro=postProcessMappedRead(inPath, anntxdb, FuSeq.SR, FuSeq.MR, FuSeq.params)
	FuSeq.SR.postPro=postProcessSplitRead(inPath, anntxdb, FuSeq.SR, FuSeq.MR, txFastaFile, FuSeq.params)
	
	myFusionFinal.MR=FuSeq.MR.postPro$myFusionFinal
	myFusionFinal.SR=FuSeq.SR.postPro$myFusionFinal

	fragmentInfo=FuSeq.MR$fragmentInfo
	FuSeq.integration=integrateFusion(myFusionFinal.MR, myFusionFinal.SR, FuSeq.params, fragmentInfo=fragmentInfo, paralog.fc.thres=2.0)
	myFusionFinal=FuSeq.integration$myFusionFinal
	
	if (nrow(myFusionFinal)==0){
		cat("\n No final fusion-gene candidates existing. If keepRData=TRUE in the parameter setting, it might be interesting to further investigate other SR and MR fusion candidates in FuSeq_process.RData.")
		myFusionExport= "# No final fusion-gene candidates existing. If keepRData=TRUE in the parameter setting, it might be interesting to further investigate other SR and MR fusion candidates in FuSeq_process.RData."
		write.table(myFusionExport, file=paste(outputDir,"/fusions.FuSeq",sep=""), col.names=FALSE, row.names = FALSE,quote = FALSE, sep="\t")
	  	save(FuSeq.params, inPath, FuSeq.integration, file=paste(outputDir,"/FuSeq_logs.RData",sep=""))
	  	#keep all RData
		if (FuSeq.params$keepRData){
		  cat("\n Saving all data of FuSeq process...")
		  save(inPath,outputDir, myFusionFinal,myFusionExport,FuSeq.params, FuSeq.MR, FuSeq.SR, FuSeq.MR.postPro, FuSeq.SR.postPro,FuSeq.integration, file=paste(outputDir,"/FuSeq_process.RData",sep=""))
		}
	}  else{
	  	myFusionOut=myFusionFinal
		myFusionOut=myFusionOut[order(myFusionOut$score, decreasing = TRUE),]
		
		selectedColumns=c("gene5","chrom5p","strand5p","brpos5.start","exonBound5p","gene3","chrom3p","strand3p","brpos3.start","exonBound3p","fusionName","SR","MR","supportRead","score")

		myFusionExport=myFusionOut[,selectedColumns]
		#rename the columns: brpos5.start --> cds.brpos5.start (break point in coding region) and exonBound5p --> brpos5.start (breaking point in exon boundary), similarly to the 3 prime side
		colnames(myFusionExport)=c("gene5","chrom5","strand5","cds.brpos5.start","brpos5","gene3","chrom3","strand3","cds.brpos3.start","brpos3","fusionName","SR.passed","MR.passed","supportRead","score")
		#get symbol names of genes
		if (exists("hgncName")){
			myFusionExport$symbol5=hgncName$hgnc_symbol[match(myFusionExport$gene5,hgncName$ensembl_gene_id)]
			myFusionExport$symbol3=hgncName$hgnc_symbol[match(myFusionExport$gene3,hgncName$ensembl_gene_id)]    
		}else{
			myFusionExport$symbol5=""
			myFusionExport$symbol3=""
		}
		#reorder the column names
		myFusionExport=myFusionExport[,c("gene5","chrom5","strand5","brpos5","cds.brpos5.start","gene3","chrom3","strand3","brpos3","cds.brpos3.start","fusionName","symbol5","symbol3","SR.passed","MR.passed","supportRead","score")]
		
		#Detect extra information here
		myFusionExport$info=rep("",nrow(myFusionExport))
		myID=unique(c(which(!is.na(myFusionOut$mitoTrans5)),which(!is.na(myFusionOut$mitoTrans3))))
		myFusionExport$info[myID]=paste(myFusionExport$info[myID],"mitochondrial translation, ",sep="")
		myID=unique(c(which(!is.na(myFusionOut$ribSub5)),which(!is.na(myFusionOut$ribSub3))))
		myFusionExport$info[myID]=paste(myFusionExport$info[myID],"cytosolic ribosomal subunit, ",sep="")
		myID=unique(c(which(!is.na(myFusionOut$ribonupro5)),which(!is.na(myFusionOut$ribonupro3))))
		myFusionExport$info[myID]=paste(myFusionExport$info[myID],"ribonucleoprotein, ",sep="")

		write.table(myFusionExport, file=paste(outputDir,"/fusions.FuSeq",sep=""), col.names=TRUE, row.names = FALSE,quote = FALSE, sep="\t")
		save(FuSeq.params, inPath, FuSeq.integration, file=paste(outputDir,"/FuSeq_logs.RData",sep=""))
		#####
		#keep all RData
		if (FuSeq.params$keepRData){
		  cat("\n Saving all data of FuSeq process...")
		  save(inPath,outputDir, myFusionFinal,myFusionExport,FuSeq.params, FuSeq.MR, FuSeq.SR, FuSeq.MR.postPro, FuSeq.SR.postPro,FuSeq.integration, file=paste(outputDir,"/FuSeq_process.RData",sep=""))
		}
		
		##### Export fasta sequence
		if (FuSeq.params$exportFasta){
			cat("\n Export supporing read sequences to files...")
			fastaOut=paste(outputDir,"/FuSeq_",sep="")
			MRinfo=getMRinfo(fusionName=as.character(myFusionFinal$fusionName),inPath=inPath,feq=FuSeq.MR$feqInfo$feq, feqFgeMap=FuSeq.MR$feqInfo$feqFgeMap, anntxdb=anntxdb,readStrands=FuSeq.params$readStrands)
			exportMappedFusionReads(inPath, readStrands=FuSeq.params$readStrands, fastaOut=fastaOut, junctInfo=MRinfo$junctInfo, fusionName=as.character(myFusionFinal$fusionName),fsizeLadder=FuSeq.MR$junctBr$fsizeLadder)
				exportSplitFusionReads(inPath, readStrands=FuSeq.params$readStrands, fastaOut=fastaOut, splitReads=FuSeq.SR$fusionGene, fusionName=as.character(myFusionFinal$fusionName))
		}		
		
	}

	cat("\n Done! \n")
}


