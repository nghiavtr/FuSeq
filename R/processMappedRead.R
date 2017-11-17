############################################################
##### process mapped reads
############################################################

processMappedRead <-function(inPath,geneAnno, anntxdb, geeqMap, FuSeq.params,feqInfo=NULL){

	cat("\n ------------------------------------------------------------------")
	cat("\n Processing mapped reads (MR) from dataset: ",inPath, " read strands:", FuSeq.params$readStrands)
	if (FuSeq.params$keepRData){
	  if (is.null(FuSeq.params$outputDir)){ 
	    FuSeq.params$outputDir=inPath
	    cat("\n No output is set, data will be saved to ", inPath)
	  }
	}
	
	if (is.null(feqInfo)){
	  feqInfo=processFEQ(inPath,geneAnno,anntxdb,geeqMap,readStrands=FuSeq.params$readStrands,chromRef=FuSeq.params$chromRef)
	  if (FuSeq.params$keepRData){
	    cat("\n Saving MR fusion candidates ...")
	    save(feqInfo,file=paste(FuSeq.params$outputDir,"/FuSeq_MR_feqInfo.RData",sep=""))
	  }
	}

	
	if (is.null(feqInfo)) return(NULL)

	#########################
	fragmentInfo=read.csv(paste(inPath,"/fragmentInfo.txt",sep=""), header =TRUE, sep="\t")
	fragDist = read.table(paste(inPath,"/fragmentDist.txt",sep=""))
	
	
	maxSharedCount=FuSeq.params$maxSharedCount;
	minGeneDist=FuSeq.params$minGeneDist;
	maxInvertedFusionCount=FuSeq.params$maxInvertedFusionCount;
	maxMRfusionFc=FuSeq.params$maxMRfusionFc;
	maxMRfusionNum=FuSeq.params$maxMRfusionNum;
	sgtMRcount=FuSeq.params$sgtMRcount;
	readStrands=FuSeq.params$readStrands
	countPropLowBound=0 # not use this one anymore 
	
	myFusionFinal=feqInfo$fgeList
	cat("\n The total number of fge candidates: ",nrow(myFusionFinal))
	
	if (readStrands=="RF" || readStrands=="UN" || readStrands=="RR"){
	  myFusionFinal$gene5p=myFusionFinal$gene2
	  myFusionFinal$gene3p=myFusionFinal$gene1
	  myFusionFinal$chrom5p=myFusionFinal$chrom2
	  myFusionFinal$chrom3p=myFusionFinal$chrom1
	  myFusionFinal$strand5p=myFusionFinal$strand2
	  myFusionFinal$strand3p=myFusionFinal$strand1
	  
	}else{
	  myFusionFinal$gene5p=myFusionFinal$gene1
	  myFusionFinal$gene3p=myFusionFinal$gene2
	  myFusionFinal$chrom5p=myFusionFinal$chrom1
	  myFusionFinal$chrom3p=myFusionFinal$chrom2
	  myFusionFinal$strand5p=myFusionFinal$strand1
	  myFusionFinal$strand3p=myFusionFinal$strand2
	}
	
	myFusionFinal$fusionName=paste(myFusionFinal$gene5p,myFusionFinal$gene3p,sep="-")

	cat("\n Start filtering ...")
	##### do some strong filters here
	if (maxSharedCount>0)	myFusionFinal=myFusionFinal[abs(myFusionFinal$supportCount-myFusionFinal$correctedCount)/myFusionFinal$supportCount <= maxSharedCount,]
	
	#remove name21 count
	if (maxInvertedFusionCount>0)	myFusionFinal=myFusionFinal[myFusionFinal$name21Count-1 <= myFusionFinal$supportCount*maxInvertedFusionCount,]
	#remove small gene distance
	if (minGeneDist > 0){
	  rmID=which(as.character(myFusionFinal$chrom1)==as.character(myFusionFinal$chrom2) & myFusionFinal$geneDist <= minGeneDist)
	  if (length(rmID)>0)
	    myFusionFinal=myFusionFinal[-rmID,]
	}
	myFusionFinal=myFusionFinal[order(myFusionFinal$supportCount, decreasing = TRUE),]
	
	
	### add gene types, keep only protein_coding genes later
	dim(geneAnno)
	matchID=match(myFusionFinal$gene1,geneAnno[,6])
	res=geneAnno[matchID,]
	colnames(res)=paste(colnames(res),"1",sep="")
	myFusionFinal=cbind(myFusionFinal,res[,c(2,4)])
	matchID=match(myFusionFinal$gene2,geneAnno[,6])
	res=geneAnno[matchID,]
	colnames(res)=paste(colnames(res),"2",sep="")
	myFusionFinal=cbind(myFusionFinal,res[,c(2,4)])
	
	#fiter low counts with multiple duplicated genes
	res=computeDupGene(myFusionFinal,dupGene.thres=-1)
	colnames(res)=c("dupGene1_f1","dupGene2_f1")
	myFusionFinal=cbind(myFusionFinal,res)
	
	rmID=(myFusionFinal$supportCount<=1) & (myFusionFinal$dupGene1_f1 > 1 | myFusionFinal$dupGene2_f1 > 1)
	myFusionFinal=myFusionFinal[!rmID,]
	
	
	#filter by dupGenes: keep maxMRfusionNum duplicated genes
	keepID1=unlist(lapply(unique(as.character(myFusionFinal$gene1)), function(g){
	  keepID=which(myFusionFinal$gene1==g)
	  if (sum(myFusionFinal$supportCount[keepID]>sgtMRcount)>maxMRfusionNum) return (NULL)
	  if (length(keepID) < maxMRfusionNum) return(keepID)
	  if(length(keepID) > maxMRfusionNum )
	    if (myFusionFinal$supportCount[keepID[maxMRfusionNum+1]]*maxMRfusionFc >= myFusionFinal$supportCount[keepID[1]]) return(NULL)
	  keepID=keepID[1:maxMRfusionNum]
	  return(keepID)
	  
	}))
	keepID2=unlist(lapply(unique(as.character(myFusionFinal$gene2)), function(g){
	  keepID=which(myFusionFinal$gene2==g)
	  if (sum(myFusionFinal$supportCount[keepID]>sgtMRcount)>maxMRfusionNum) return (NULL)
	  if (length(keepID) < maxMRfusionNum) return(keepID)
	  if(length(keepID) > maxMRfusionNum )
	    if (myFusionFinal$supportCount[keepID[maxMRfusionNum+1]]*maxMRfusionFc >= myFusionFinal$supportCount[keepID[1]]) return(NULL)
	  keepID=keepID[1:maxMRfusionNum]
	  return(keepID)
	}))
	#  keepID=unique(c(keepID1,keepID2))
	keepID=intersect(keepID1,keepID2)
	myFusionFinal=myFusionFinal[keepID,]
	myFusionFinal=myFusionFinal[order(myFusionFinal$supportCount, decreasing = TRUE),]
	
	
	cat("\n The number of remaining fge candidates: ",nrow(myFusionFinal))
	######### compute Count of fge from ftx
	cat("\n Get the sum count of ftx")
	
	feqRaw1=feqInfo$feqRaw[feqInfo$feqRaw$Read==1,]
	feqRaw2=feqInfo$feqRaw[feqInfo$feqRaw$Read==2,]
	
	myfeqID=lapply(as.character(myFusionFinal$name12), function(mykey) feqInfo$feqFgeMap[[mykey]])
	myfeqIDSize=sapply(myfeqID, length)
	myfeqIDName=unlist(apply(cbind(as.character(myFusionFinal$name12),myfeqIDSize),1, function(x) rep(x[1],x[2])))
	myfeqID=unlist(myfeqID)
	
	ufeqID=unique(myfeqID)
	mydat=feqRaw1[!is.na(match(feqRaw1$Feq,ufeqID)),]
	tx1Num=tapply(mydat[,3],mydat[,4],sum)
	mydat=feqRaw2[!is.na(match(feqRaw2$Feq,ufeqID)),]
	tx2Num=tapply(mydat[,3],mydat[,4],sum)/2
	mytxCount=tx1Num*tx2Num*feqInfo$feq[as.integer(names(tx2Num))]
	matchID=match(myfeqID,as.integer(names(tx2Num)))
	
	myres=tapply(mytxCount[matchID],myfeqIDName,sum)
	myres=myres[match(as.character(myFusionFinal$name12),names(myres))]
	
	#get the aggregation counts from tx
	myFusionFinal$aggtxCount=myres
	myFusionFinal$countProp=myFusionFinal$aggtxCount/myFusionFinal$supportCount
	
	#get the number of tx of gene1
	myres=tapply(tx1Num[matchID],myfeqIDName,sum)
	myres=myres[match(as.character(myFusionFinal$name12),names(myres))]
	myFusionFinal$tx1Num=myres
	#get the number of tx of gene2
	myres=tapply(tx2Num[matchID],myfeqIDName,sum)
	myres=myres[match(as.character(myFusionFinal$name12),names(myres))]
	myFusionFinal$tx2Num=myres
	
	
	cat("\n Continue the filtering ...")
	
	#filter by count
	myFusionFinal=myFusionFinal[myFusionFinal$supportCount >=FuSeq.params$minMR,]
	
	#filter by countProp
	myFusionFinal=myFusionFinal[myFusionFinal$countProp>countPropLowBound,] 
	
	cat("\n The number of remaining fge candidates: ",nrow(myFusionFinal))
	
	#do biological filters
	cat("\n Filter by biological features... ")
	bioFilter.res=doBiologicalFilter(myFusionFinal, chromRef=FuSeq.params$chromRef, onlyProteinCodingGenes=FuSeq.params$onlyProteinCodingGenes, doFilter=TRUE)
	myFusionFinal=bioFilter.res
	
	cat("\n The number of remaining fge candidates: ",nrow(myFusionFinal))
	#sequence similarity between two genes
	cat("\n Filter by sequence similarity... ")
	seqHmlog=NULL
	for (i in 1:nrow(myFusionFinal)){
	  res=intersect(geeqMap[[as.character(myFusionFinal$gene2[i])]],geeqMap[[as.character(myFusionFinal$gene1[i])]])
	  if (length(res) > 0) seqHmlog=c(seqHmlog,1) else (seqHmlog=c(seqHmlog,0))
	}
	myFusionFinal$seqHmlog=seqHmlog
	
	#do filter
	myFusionFinal=myFusionFinal[myFusionFinal$seqHmlog==0,]
	cat("\n The number of remaining fge candidates: ",nrow(myFusionFinal))
	
	
	#detect junction breaks
	cat("\n Detect junction breaks... ")
	junctBr=detectJunctionBreaks(myFusionFinal,inPath, feqInfo$feq,feqInfo$feqFgeMap, anntxdb, readStrands=FuSeq.params$readStrands)
	myFusionFinal=junctBr$myFusionFinal
	
	myFusionFinal=myFusionFinal[order(myFusionFinal$nondupCount, decreasing = TRUE),]
	#filter by nondupCount
	myFusionFinal=myFusionFinal[myFusionFinal$nondupCount>=FuSeq.params$minNonDupMR,]
	
	#filter by junction distance
	myFusionFinal=myFusionFinal[myFusionFinal$junctDist>FuSeq.params$minJunctionDist,]
	
	if (nrow(myFusionFinal) == 0){ 
	  cat("\n The number of final fge candidates: ",nrow(myFusionFinal))
	  return(NULL)
	}

	#get coverage
	myFusionFinal$cover5=myFusionFinal$genebrpos5.rg/myFusionFinal$flen5
	myFusionFinal$cover3=myFusionFinal$genebrpos3.rg/myFusionFinal$flen3
	

	cat("\n The number of final fge candidates: ",nrow(myFusionFinal))

	res=list(myFusionFinal=myFusionFinal,junctBr=junctBr,feqInfo=feqInfo,fragmentInfo=fragmentInfo,fragDist=fragDist)
	return(res)
}
