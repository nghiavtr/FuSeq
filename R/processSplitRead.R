############################################################
##### process split reads
############################################################

processSplitRead <-function(inPath,geneAnno, anntxdb, geeqMap, txFastaFile, FuSeq.params){
  cat("\n ------------------------------------------------------------------")
  cat("\n Processing split reads (SR) from dataset: ",inPath, " read strands:", FuSeq.params$readStrands)
  
  ##### get fragment information
  fragmentInfo=read.csv(paste(inPath,"/fragmentInfo.txt",sep=""), header =TRUE, sep="\t")
  fragDist = read.table(paste(inPath,"/fragmentDist.txt",sep=""))
  fragRg=fragDist[fragDist[,2]>0,1]
  flen.min=min(fragRg)
  flen.max=max(fragRg)
  ##### find split reads
  cat("\n Get split reads ...")
  splitReads=fsizes=NULL
  frfiles=list.files(inPath,paste("splitReadInfo_*",sep=""))
  for (i in 1:length(frfiles)){
    tmpDat=read.csv(paste(inPath,"/",frfiles[i],sep=""), header =FALSE, sep="\t")
    splitReads=rbind(splitReads,tmpDat)
    fsizes=c(fsizes,nrow(tmpDat))
  }
  fsizeLadder=cumsum(fsizes)
  colnames(splitReads)=c("header","readType","direction","front_tx","front_gene","front_hitpos","front_querypos","front_len","back_tx","back_gene","back_hitpos","back_querypos","back_len","matchedGene","matchedDirect","matchedPos")
  splitReads$tx12=paste(splitReads$front_tx,splitReads$back_tx,sep="-")
  splitReads$tx21=paste(splitReads$back_tx,splitReads$front_tx,sep="-")
  splitReads$name12=paste(splitReads$front_gene,splitReads$back_gene,sep="-")
  splitReads$name21=paste(splitReads$back_gene,splitReads$front_gene,sep="-")
  
  
  ############# starting process
  fusionGene=splitReads

  cat("\n Extract other biological information...")
  #add few biological information
  res=select(anntxdb, keys=as.character(fusionGene$front_gene), columns=c("GENEID","TXCHROM","TXSTRAND"), keytype = "GENEID")
  colnames(res)=c("GENEID","chrom1","strand1")
  fusionGene=cbind(fusionGene,res[,-1])
  res=select(anntxdb, keys=as.character(fusionGene$back_gene), columns=c("GENEID","TXCHROM","TXSTRAND"), keytype = "GENEID")
  colnames(res)=c("GENEID","chrom2","strand2")
  fusionGene=cbind(fusionGene,res[,-1])
  fusionGene=chromFilter(fusionGene) 
  
  fusionGene$gene1=fusionGene$front_gene
  fusionGene$gene2=fusionGene$back_gene
  
  
  matchID=match(fusionGene$gene1,geneAnno[,6])
  res=geneAnno[matchID,]
  colnames(res)=paste(colnames(res),"1",sep="")
  fusionGene=cbind(fusionGene,res[,c(2,4)])
  matchID=match(fusionGene$gene2,geneAnno[,6])
  res=geneAnno[matchID,]
  colnames(res)=paste(colnames(res),"2",sep="")
  fusionGene=cbind(fusionGene,res[,c(2,4)])
  #filter by protein coding
  keepID=which(fusionGene$geneType1=="protein_coding" & fusionGene$geneType2=="protein_coding")
  fusionGene=fusionGene[keepID,]
  
  
  ##### get supporting count
  res=table(fusionGene$name12)
  fusionGene$supportCount=res[match(fusionGene$name12,names(res))]

  ##### count of name21
  matchID=match(fusionGene$tx12, fusionGene$tx21)
  tx21Count=fusionGene$supportCount[matchID]
  tx21Count[is.na(tx21Count)]=0
  fusionGene$tx21Count=tx21Count
  
  ##### shared count
  res=table(as.character(fusionGene$header))
  length(res)
  fusionGene$readProp=1/res[match(as.character(fusionGene$header),names(res))]
  
  res=table(fusionGene$tx12)
  fusionGene$tx12Count=res[match(fusionGene$tx12,names(res))]
  adjtx12Count=tapply(as.double(fusionGene$readProp),as.character(fusionGene$tx12),sum)
  fusionGene$adjtx12Count=adjtx12Count[match(as.character(fusionGene$tx12),names(adjtx12Count))]
  
  #gene distance
  geneDist=computeGeneDistance(fusionGene,anntxdb,minGeneDist=-1)
  fusionGene$geneDist=geneDist
  
  myFusion=fusionGene
  
  #filter by gene distance
  rmID=which(as.character(myFusion$chrom1)==as.character(myFusion$chrom2) &  myFusion$geneDist <= FuSeq.params$minGeneDist)
  if (length(rmID)>0) myFusion=myFusion[-rmID,]
  
  #filter by overlapping sequences
  myFusion$brol=myFusion$front_querypos+myFusion$front_len-myFusion$back_querypos
  myFusionFP=myFusion[myFusion$brol>10,]
  myFusion=myFusion[myFusion$brol<=10,]
  
  FPlist=unique(as.character(myFusionFP$name12))
  matchID=match(as.character(myFusion$name12),FPlist)
  myFusion=myFusion[which(is.na(matchID)),]

  res=table(myFusion$name12)
  myFusion$supportCount2=res[match(myFusion$name12,names(res))]
  #transcript level
  res=table(myFusion$tx12)
  myFusion$tx12Count2=res[match(myFusion$tx12,names(res))]
  #remove all ftx not satisfying the sequence overlapping
  keepID=which(myFusion$tx12Count2==myFusion$tx12Count)
  myFusion=myFusion[keepID,]
  
  
  myFusion$flen=rep(0,nrow(myFusion))
  fwID=which(myFusion$direction==0 | myFusion$direction==3)
  rcID=which(myFusion$direction==1 | myFusion$direction==4)
  #forward
  myFusion$flen[fwID]=myFusion$back_querypos[fwID] +  myFusion$matchedPos[fwID]-myFusion$back_hitpos[fwID] + fragmentInfo$readlen
  #rc
  myFusion$flen[rcID]=(myFusion$front_hitpos[rcID] -myFusion$front_querypos[rcID]- myFusion$matchedPos[rcID] ) + fragmentInfo$readlen
  myFusion=myFusion[myFusion$flen>=flen.min,]
  myFusion=myFusion[myFusion$flen<=flen.max,]
  #mySR$front_hitpos+mySR$front_len-mySR$front_querypos
  mytest=pnorm(myFusion$flen, mean=fragmentInfo$fragLengthMean, sd=fragmentInfo$fragLengthSd)
  myFusion$flenTest=mytest
  myFusion=myFusion[myFusion$flenTest>=0.001,]
  myFusion=myFusion[myFusion$flenTest<=0.999,]
  

  ### find duplicate tx
  myDup=duplicated(myFusion$tx12)
  myFusionNoDup=myFusion[!myDup,]
  
  duptx1=table(myFusionNoDup$front_tx)
  duptx2=table(myFusionNoDup$back_tx)
  myFusion=cbind(myFusion,as.integer(duptx1[match(as.character(myFusion$front_tx),names(duptx1))]))
  myFusion=cbind(myFusion,as.integer(duptx2[match(as.character(myFusion$back_tx),names(duptx2))]))
  colnames(myFusion)[c(ncol(myFusion)-1,ncol(myFusion))]=c("duptx1","duptx2")
  
 
  ####### check canonical splicing sites
  cat("\n Remaining fusion reads: ",nrow(myFusion))
  cat("\n Check canonical splicing sites, it takes time ...")
  
  #####preparing annotation
  library(Biostrings)
  fasta = readDNAStringSet(txFastaFile)
  fasta_txnames=sapply(names(fasta), function(x) unlist(strsplit(x," "))[1])
  
  shrinkLen=5
  myFusion$front_brpos=myFusion$front_hitpos+myFusion$front_len-1
  myFusion$brchposEx5=-1
  myFusion$back_brpos=myFusion$back_hitpos-(myFusion$back_querypos-myFusion$front_querypos-myFusion$front_len)-1
  myFusion$brchposEx3=-1
  

  myFusion$GTATCEnd=myFusion$ATATCEnd=myFusion$ATEnd=myFusion$GCEnd=myFusion$GTEnd=myFusion$ssExEnd=myFusion$ssExEndGe=rep(-1,nrow(myFusion))
  myFusion$AAStart=myFusion$ATStart=myFusion$ACStart=myFusion$AGStart=myFusion$ssExStart=myFusion$ssExStartGe=rep(-1,nrow(myFusion))
  ##### check splicing
  cat("\n Checking in 5 prime site...")
  
  #to speedup we do in blocks
  blockSize=20000
  blockNum=trunc(nrow(myFusion)%/%blockSize) + ifelse(nrow(myFusion)%%blockSize>0,1,0)
  #sorted by front_tx
  myFusion=myFusion[order(myFusion$front_tx),]
  
  mytime <- system.time({
  
  for (blockID in 1:blockNum){
    cat("\n",blockID," blocks processed")
    
    block.keepID=c(((blockID-1)*blockSize+1):(blockID*blockSize))
    block.keepID=block.keepID[block.keepID<=nrow(myFusion)]
    myFusionBlock=myFusion[block.keepID,]
    
    exonInfo=select(anntxdb, keys=unique(c(as.character(myFusionBlock$front_gene))), columns=c("TXNAME","EXONID","EXONSTART","EXONEND","TXSTRAND"), keytype = "GENEID")
      frontTxList=as.character(myFusionBlock$front_tx)
      frontTxSet=unique(frontTxList)
      cat("\n Total transcripts ",length(frontTxSet))
      for (i in 1:length(frontTxSet)){
        
        keepID=which(frontTxList==frontTxSet[i])
        mySR=myFusionBlock[keepID,]
        txname=frontTxSet[i]
        
        matchID=match(txname,fasta_txnames)
        mytxFasta=fasta[matchID]
        txExonMat=exonInfo[exonInfo$TXNAME==txname,]
        
        if (txExonMat$TXSTRAND[1]=="+") {
          txExonMat=txExonMat[order(txExonMat$EXONSTART),]#sort exons by increasing order for forward strand
          exonlen=txExonMat$EXONEND-txExonMat$EXONSTART+1
          exonlenCumSum=cumsum(exonlen)
          txExonMat$EXONTXSTART=c(0,exonlenCumSum[-length(exonlenCumSum)])
          txExonMat$EXONTXEND=exonlenCumSum-1
        }else{
          txExonMat=txExonMat[order(txExonMat$EXONEND, decreasing=TRUE),]#sort exons by decreasing order for reverse strand
          exonlen=txExonMat$EXONEND-txExonMat$EXONSTART+1
          exonlenCumSum=cumsum(exonlen)
          txExonMat$EXONTXSTART=c(0,exonlenCumSum[-length(exonlenCumSum)])
          txExonMat$EXONTXEND=exonlenCumSum-1
        }
        
        myDiff=NULL
        myDonor5=myDonor=NULL
        if (txExonMat$TXSTRAND[1]=="-") diff=unlist(lapply(mySR$front_brpos,function(x) min(abs(x-txExonMat$EXONTXSTART)))) else diff=unlist(lapply(mySR$front_brpos,function(x) min(abs(x-txExonMat$EXONTXEND))))
        myDiff=diff
        
        for (k in 1:shrinkLen){
          mybrpos=mySR$front_brpos-k+1
          donorSeq=substring(mytxFasta,mybrpos-1,mybrpos)
          myDonor=cbind(myDonor,donorSeq)
          
          donorSeq5=substring(mytxFasta,mybrpos-4,mybrpos)
          myDonor5=cbind(myDonor5,donorSeq5)
          
        }
        exEnd=myDiff
        myFusionBlock$ssExEnd[keepID]=exEnd
        myDonorCk=apply(myDonor,1,function(x) max(x=="GT"))
        myFusionBlock$GTEnd[keepID]=myDonorCk
        myDonorCk=apply(myDonor,1,function(x) max(x=="GC"))
        myFusionBlock$GCEnd[keepID]=myDonorCk
        myDonorCk=apply(myDonor,1,function(x) max(x=="AT"))
        myFusionBlock$ATEnd[keepID]=myDonorCk
        
        myDonorCk=apply(myDonor5,1,function(x) max(x=="ATATC"))
        myFusionBlock$ATATCEnd[keepID]=myDonorCk

        myDonorCk=apply(myDonor5,1,function(x) max(x=="GTATC"))
        myFusionBlock$GTATCEnd[keepID]=myDonorCk
        
        
        myFusionBlock$brchposEx5[keepID]=convertChrPos(txname=txname,txpos=myFusionBlock$front_brpos[keepID],txExonMat=txExonMat)      
        #check exon at gene-level
        geExonMat=exonInfo[exonInfo$GENEID==mySR$gene1[1],]
        myDiff=NULL
        #diff=unlist(lapply(myFusionBlock$brchposEx5[keepID],function(x) min(abs(x-c(geExonMat$EXONSTART,geExonMat$EXONEND))))) #this one is less strict where we allow a shift between exons. Results of this will be concordant with the results from tx checking
        if (geExonMat$TXSTRAND[1]=="-") diff=unlist(lapply(myFusionBlock$brchposEx5[keepID],function(x) min(abs(x-geExonMat$EXONSTART)))) else diff=unlist(lapply(myFusionBlock$brchposEx5[keepID],function(x) min(abs(x-geExonMat$EXONEND)))) #this one is more strict, so the result might be not concordant with the tx checking
        myDiff=diff
        
        
        myFusionBlock$ssExEndGe[keepID]=myDiff
        
        if (i %% 1000 ==0) cat("\n",i," transcripts processed")
      }
    #update results
    myFusion[block.keepID,]=myFusionBlock
  }

  })
    mytime
    
    
  
  cat("\n\n Checking in 3 prime site...")
  
  #to speedup we do in blocks
  blockSize=20000
  blockNum=trunc(nrow(myFusion)%/%blockSize) + ifelse(nrow(myFusion)%%blockSize>0,1,0)
  #sorted by back_tx
  myFusion=myFusion[order(myFusion$back_tx),]
  
  mytime <- system.time({
    
    for (blockID in 1:blockNum){
      cat("\n",blockID," blocks processed")
      
      block.keepID=c(((blockID-1)*blockSize+1):(blockID*blockSize))
      block.keepID=block.keepID[block.keepID<=nrow(myFusion)]
      myFusionBlock=myFusion[block.keepID,]
      
      exonInfo=select(anntxdb, keys=unique(c(as.character(myFusionBlock$back_gene))), columns=c("TXNAME","EXONID","EXONSTART","EXONEND","TXSTRAND"), keytype = "GENEID")
      
        backTxList=as.character(myFusionBlock$back_tx)
        backTxSet=unique(backTxList)
        cat("\n Total transcripts ",length(backTxSet))
        for (i in 1:length(backTxSet)){
          keepID=which(backTxList==backTxSet[i])
          mySR=myFusionBlock[keepID,]
          txname=backTxSet[i]
          
          matchID=match(txname,fasta_txnames)
          mytxFasta=fasta[matchID]
          txExonMat=exonInfo[exonInfo$TXNAME==txname,]
          if (txExonMat$TXSTRAND[1]=="+") {
            txExonMat=txExonMat[order(txExonMat$EXONSTART),]#sort exons by increasing order for forward strand
            exonlen=txExonMat$EXONEND-txExonMat$EXONSTART+1
            exonlenCumSum=cumsum(exonlen)
            txExonMat$EXONTXSTART=c(0,exonlenCumSum[-length(exonlenCumSum)])
            txExonMat$EXONTXEND=exonlenCumSum-1
          }else{
            txExonMat=txExonMat[order(txExonMat$EXONEND, decreasing=TRUE),]#sort exons by decreasing order for reverse strand
            exonlen=txExonMat$EXONEND-txExonMat$EXONSTART+1
            exonlenCumSum=cumsum(exonlen)
            txExonMat$EXONTXSTART=c(0,exonlenCumSum[-length(exonlenCumSum)])
            txExonMat$EXONTXEND=exonlenCumSum-1
          }
          
          myDiff=NULL
          myAcceptor=NULL
          if (txExonMat$TXSTRAND[1]=="+") diff=unlist(lapply(mySR$back_brpos,function(x) min(abs(x-txExonMat$EXONTXSTART)))) else diff=unlist(lapply(mySR$back_brpos,function(x) min(abs(x-txExonMat$EXONTXEND))))
          myDiff=diff
          
          for (k in 1:shrinkLen){
            mybrpos=mySR$back_brpos+k
            acceptorSeq=substring(mytxFasta,mybrpos-1,mybrpos)
            myAcceptor=cbind(myAcceptor,acceptorSeq)
          }
          #exStart=rowMin(myDiff)
          exStart=myDiff
          myFusionBlock$ssExStart[keepID]=exStart

          myAcceptorCk=apply(myAcceptor,1,function(x) max(x=="AG"))
          myFusionBlock$AGStart[keepID]=myAcceptorCk
          myAcceptorCk=apply(myAcceptor,1,function(x) max(x=="AC"))
          myFusionBlock$ACStart[keepID]=myAcceptorCk
          myAcceptorCk=apply(myAcceptor,1,function(x) max(x=="AT"))
          myFusionBlock$ATStart[keepID]=myAcceptorCk
          myAcceptorCk=apply(myAcceptor,1,function(x) max(x=="AA"))
          myFusionBlock$AAStart[keepID]=myAcceptorCk
          
          
          myFusionBlock$brchposEx3[keepID]=convertChrPos(txname=txname,txpos=myFusionBlock$back_brpos[keepID],txExonMat=txExonMat)     
          #check exon at gene-level
          geExonMat=exonInfo[exonInfo$GENEID==mySR$gene2[1],]
          myDiff=NULL
          #diff=unlist(lapply(myFusionBlock$brchposEx3[keepID],function(x) min(abs(x-c(geExonMat$EXONSTART,geExonMat$EXONEND))))) #this one is less strict where we allow a shift between exons. Results of this will be concordant with the results from tx checking
          if (geExonMat$TXSTRAND[1]=="+") diff=unlist(lapply(myFusionBlock$brchposEx3[keepID],function(x) min(abs(x-geExonMat$EXONSTART)))) else diff=unlist(lapply(myFusionBlock$brchposEx3[keepID],function(x) min(abs(x-geExonMat$EXONEND)))) #this one is more strict, so the result might be not concordant with the tx checking
          myDiff=diff
          
          myFusionBlock$ssExStartGe[keepID]=myDiff
          
          if (i %% 1000 ==0) cat("\n",i," processed")
        }
        
 
      
      #update results
      myFusion[block.keepID,]=myFusionBlock
    }
    
  })
  mytime
  
  myFusion$junctDist=abs(myFusion$brchposEx5-myFusion$brchposEx3)
  
  myFusion$ssStart=myFusion$ssEnd=rep(-1,nrow(myFusion))
  myFusion$ssEnd=ifelse(abs(myFusion$GCEnd)>0,4,myFusion$ssEnd)
  myFusion$ssEnd=ifelse(abs(myFusion$GTEnd)>0,3,myFusion$ssEnd)
  #myFusion$ssEnd=ifelse(abs(myFusion$ssCkEnd)>0,3,myFusion$ssEnd)
  myFusion$ssEnd=ifelse(abs(myFusion$ssExEndGe)<shrinkLen,2,myFusion$ssEnd)
  myFusion$ssEnd=ifelse(abs(myFusion$ssExEnd)<shrinkLen,1,myFusion$ssEnd)
  
  myFusion$ssStart=ifelse(abs(myFusion$AGStart)>0,3,myFusion$ssStart)
  #myFusion$ssStart=ifelse(abs(myFusion$ssCkStart)>0,3,myFusion$ssStart)
  myFusion$ssStart=ifelse(abs(myFusion$ssExStartGe)<shrinkLen,2,myFusion$ssStart)
  myFusion$ssStart=ifelse(abs(myFusion$ssExStart)<shrinkLen,1,myFusion$ssStart)
  
  
  
  #keep only split reads passed the splicing site condition
  myFusionTmp=myFusion
  myFusionTmp=myFusionTmp[myFusionTmp$ssStart>0,]
  myFusionTmp=myFusionTmp[myFusionTmp$ssEnd>0,]

  rmID=which(myFusionTmp$ssEnd>=3 & myFusionTmp$ssStart>=3)
  if (length(rmID)>0)  myFusionTmp=myFusionTmp[-rmID,]
  
   cat("\n\n Checking possible paralogs...")

  #remove paralogs from database
  rmID=unlist(lapply(c(1:nrow(myFusionTmp)), function(x){
    par1=c(as.character(myFusionTmp$gene1[x]),geneParalog[which(geneParalog[,1]==as.character(myFusionTmp$gene1[x])),2])
    par2=c(as.character(myFusionTmp$gene2[x]),geneParalog[which(geneParalog[,1]==as.character(myFusionTmp$gene2[x])),2])
    if (length(par1)>0 & length(par2)>0) return(length(intersect(par1,par2))>0)
    return(FALSE)
  }))
  
  myFusionTmp$paralog=rep(0,nrow(myFusionTmp))
  myFusionTmp$paralog[rmID]=1
  myFusionTmp=myFusionTmp[myFusionTmp$paralog==0,]
  

  myFusionFinal=myFusionTmp
	res=list(myFusionFinal=myFusionFinal,myFusion=myFusion,fusionGene=fusionGene, splitReads=splitReads,fragmentInfo=fragmentInfo,fragDist=fragDist,myFusionFP=myFusionFP)
	return(res)
}




