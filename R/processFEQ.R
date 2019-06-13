############################################################
### 13May2019/Nghia: 
#- revise codes to improve speed
#- allow to get input from multiple feq files
#####process fusion equivalence classes to generate initial fusion gene candidates
processFEQ <-function(inPath,anntxdb,readStrands="UN",chromRef=as.character(c(1:22,"X","Y"))){
  cat("\n Read fusion equivalence classes")
  ##### read fragment information
  fragmentInfo=read.csv(paste(inPath,"/fragmentInfo.txt",sep=""), header =TRUE, sep="\t")
  libsize=fragmentInfo[1,5]
  #read fusion-equivalence classes
  #if (readStrands=="ALL") readStrands=c("RR","FF","RF","FR")  
  feqRaw=NULL
  eqNum=0
  for(rf in readStrands){
    fnIn=paste(inPath,"/feq_",rf,".txt",sep="")
    cat("\n Reading ",fnIn)
    feqR=read.csv(fnIn, header =TRUE, sep="\t",as.is=TRUE)
    maxEq=max(feqR$Feq)
    feqR=cbind(feqR,readStrands=rep(rf,nrow(feqR)))
    feqR$Feq=feqR$Feq+eqNum
    #update eqNum
    eqNum=eqNum+maxEq
    #add up
    feqRaw=rbind(feqRaw,feqR)
  }

 ### merging the same feq
#  #get unique feq
#  feqRaw_TxSet=paste(feqRaw$Read,feqRaw$Transcript,sep=":")
#  eqClass = tapply(feqRaw_TxSet,feqRaw$Feq, function(x){paste(sort(x),collapse='::')})  
#  myDup=duplicated(feqRaw$Feq)  
#  eqCount=feqRaw$Count[!myDup]
#  names(eqCount)=feqRaw$Feq[!myDup]
#  eqCount2=tapply(eqCount,eqClass,sum)
#  feqRaw_eqClass=eqClass[match(feqRaw$Feq,as.integer(names(eqClass)))]
#  feqRaw$Count_raw=feqRaw$Count
#  feqRaw$Count=eqCount2[match(feqRaw_eqClass,names(eqCount2))]
#  #update new index
#  feqRaw$Feq_raw=feqRaw$Feq
#  newFeq=c(1:length(eqCount2))
#  feqRaw$Feq=newFeq[match(feqRaw_eqClass,names(eqCount2))]
#  #reindex again
#  newFeqUni=unique(feqRaw$Feq)
#  feqRaw$Feq=match(feqRaw$Feq,newFeqUni)  
#  feqRaw_eqClass=NULL
#  feqRaw_TxSet=NULL
  #### back to original FuSeq
  #generate feq  
  myDup=duplicated(feqRaw$Feq)  
  feq=feqRaw$Count[!myDup]
  names(feq)=feqRaw$Feq[!myDup]
  cat("\n The total number of fusion equivalence classes: ",length(feq))  
  if (length(feq) == 0) return(NULL)
  
  cat("\n Create the maps between feq and fge")
  #get a map between tx and gene
  txToGene=select(anntxdb, keys=as.character(feqRaw[,1]), columns=c("GENEID","TXCHROM"), keytype = "TXNAME")
  feqRaw$GENEID=txToGene$GENEID
  feqGene=feqRaw[,c("Read","Feq","GENEID")]
  #remove duplicates
  myDup=duplicated(feqGene)
  feqGene=feqGene[!myDup,]
  feqGene=feqGene[order(feqGene$Feq),] #they must be in increasing order 1,2...
  feqGene1=feqGene[feqGene$Read==1,]
  feqGene2=feqGene[feqGene$Read==2,]
  #remove duplicates of txToGene
  myDup=duplicated(txToGene)
  txToGene=txToGene[!myDup,]
  ##### do indexing feq and fge
  feqGene1Num=table(feqGene1$Feq)
  feqGene2Num=table(feqGene2$Feq)

 # NOTE: switch the number of tx between read1 and read2 to get all combinations
  feqGene1$Read=feqGene2Num[feqGene1$Feq]
  feqGene2$Read=feqGene1Num[feqGene2$Feq]

  gene1list=rep(feqGene1[,3],feqGene1[,1])
  res=tapply(feqGene2$GENEID,feqGene2$Feq,c)
  res=data.frame(y=feqGene1Num,x=res)
  gene2list=rep(res[,3],res[,2])
  gene2list=unlist(gene2list)
  feqID=rep(seq_along(feqGene1Num),feqGene1Num*feqGene2Num)  
  fgeneNames=paste(gene1list,gene2list,sep="-")
  
  fgeneNames_uni=unique(fgeneNames)
  fgeneNames_uni_ID=seq(length(fgeneNames_uni))
  fgeneNames_ID=fgeneNames_uni_ID[match(fgeneNames,fgeneNames_uni)]

  feqFgeMap=tapply(feqID,fgeneNames_ID,c)
  feqFgeMap=lapply(feqFgeMap,function(x) as.integer(unlist(x)))
  names(feqFgeMap)=fgeneNames_uni[match(as.integer(names(feqFgeMap)),fgeneNames_uni_ID)]

  #length(feqFgeMap) # should be equal to the number of fge
  #create fge
  fge=data.frame(fgeID=c(1:length(feqFgeMap)), name12=names(feqFgeMap),stringsAsFactors = FALSE)
   matchID=match(fgeneNames,fge$name12)
  fgeFeqMap=tapply(matchID,feqID,c)
  fgeFeqMap=lapply(fgeFeqMap,function(x) as.integer(unlist(x)))
  #length(fgeFeqMap) #shoud be equal to the number of feq 
  
  ##### add features to fge
#  res=lapply(fge$name12, function(x) unlist(strsplit(x,"-")))
  res=strsplit(fge$name12,"-")
  res=do.call(rbind,res)
  colnames(res)=c("gene1","gene2")
  fge=cbind(fge,res)
  fge$name21=paste(fge$gene2,fge$gene1,sep="-")
  
  ###########
  cat("\n Get the number of supporting reads")
  ### find raw counts of fusion genes
  fusionGene=fge
  #sum-up all counts of feq
  supportCount=rep(0,nrow(fge))
  for (feqID in 1:length(feq)){
    feqCount=feq[feqID]
    gid=unlist(fgeFeqMap[feqID])
    supportCount[gid]=supportCount[gid]+feqCount
  }
  fusionGene$supportCount=supportCount

  ####
 cat("\n Get the number of uniquely mapped supporting reads")
 feqSize=sapply(fgeFeqMap,length)
 uniFeq=which(feqSize==1)
 uniCount=rep(0,nrow(fge))
#  for (feqID in 1:length(feq)){
   for (feqID in uniFeq){
    feqCount=feq[feqID]
    gid=unlist(fgeFeqMap[feqID])
    uniCount[gid]=uniCount[gid]+feqCount
  }
  fusionGene$uniCount=uniCount  

  ### estimate counts for fusion genes
  cat("\n Correct the number of supporting reads")
  alp0=rep(sum(feq)/nrow(fusionGene),nrow(fusionGene))
  alpOut=rep(0,nrow(fusionGene))
  
  res=estimateCountEM(alp0,alpOut,feq,fgeFeqMap,itNum=1,alpDiff.thres=0.01)
  
  fusionGene$correctedCount=res$correctedCount
  #sort by corrected count
  fusionGene=fusionGene[order(fusionGene$correctedCount,decreasing=TRUE),] 
  
  #######
  cat("\n Keep only candidates in the selected chromosomes")
  
  cat("\n Extract other biological information...")
  #add few biological information
  allgenes=genes(anntxdb,single.strand.genes.only=FALSE)
  allgenes.info=select(anntxdb, keys=names(allgenes), columns=c("TXCHROM","TXSTRAND"), keytype = "GENEID")
  ##Check: it is sure that for GRCh37.75, none of geneID from 1-22,X,Y overlapping with the rest
  #allgenes.info1=allgenes.info[allgenes.info$TXCHROM %in% chromRef,]
  #allgenes.info2=allgenes.info[!allgenes.info$TXCHROM %in% chromRef,]
  #sum(!is.na(match(allgenes.info2$GENEID, allgenes.info1$GENEID)))
  
  #keep only candidates from selected chromosomes in chromRef
  allgenes.info=allgenes.info[allgenes.info$TXCHROM %in% chromRef,]
  
  #res=select(anntxdb, keys=as.character(fusionGene$gene1), columns=c("GENEID","TXCHROM","TXSTRAND"), keytype = "GENEID")
  #ptm <- proc.time()
  res=allgenes.info[match(as.character(fusionGene$gene1),allgenes.info$GENEID),c("GENEID","TXCHROM","TXSTRAND")]
  #proc.time() - ptm
  #select(anntxdb, keys=as.character(fusionGene$gene1), columns=c("GENEID","TXCHROM","TXSTRAND"), keytype = "GENEID")
  colnames(res)=c("GENEID","chrom1","strand1")
  fusionGene=cbind(fusionGene,res[,-1])
  #res=select(anntxdb, keys=as.character(fusionGene$gene2), columns=c("GENEID","TXCHROM","TXSTRAND"), keytype = "GENEID")
  res=allgenes.info[match(as.character(fusionGene$gene2),allgenes.info$GENEID),c("GENEID","TXCHROM","TXSTRAND")]
  colnames(res)=c("GENEID","chrom2","strand2")
  fusionGene=cbind(fusionGene,res[,-1])
  
  #do filter by chromosomes
  fusionGene=chromFilter(fusionGene,chromRef=chromRef) 
  
  
  ##### generate features to detect fusion genes
  cat("\n Extract extra features...")
  #Fusion-genes sharing the same set of fusion equivalence classes
  matchID=match(as.character(fusionGene$name12),names(feqFgeMap))
  feqDup=duplicated(feqFgeMap[matchID], fromLast=FALSE) + duplicated(feqFgeMap[matchID], fromLast=TRUE)
  fusionGene$feqDup=feqDup
  
  #compute cpm
  fusionGene$cpm=fusionGene$supportCount*10^6/libsize
  #distance between two genes >= minGeneDist 
  geneDist=computeGeneDistance(fusionGene,anntxdb)
  fusionGene$geneDist=geneDist
  
  #get reversed fusion (3-5) count
  name21Count=computeReversedFusionCount(fusionGene)
  fusionGene$name21Count=name21Count
  
  #compute dupGene
  res=computeDupGene(fusionGene,dupGene.thres=-1)
  fusionGene=cbind(fusionGene,res)
  
  res=list(fgeList=fusionGene,fge=fge,feq=feq,fgeFeqMap=fgeFeqMap,feqFgeMap=feqFgeMap,feqRaw=feqRaw)
  return(res)
  
}
