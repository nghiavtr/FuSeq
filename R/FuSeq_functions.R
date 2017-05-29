#Date:19/05/2017
#- Fix and clean codes
###################
##### Several functions for fusion gene detection

chromFilter <- function(myfusionTx,chromRef=NULL){
  ##### Filtering by chromosomes
  if (is.null(chromRef)) chromRef=as.character(c(1:22,"X","Y"))
  keepID=which(!is.na(match(as.character(myfusionTx$chrom1),chromRef)))
  myfusionTx=myfusionTx[keepID,]
  keepID=which(!is.na(match(as.character(myfusionTx$chrom2),chromRef)))
  myfusionTx=myfusionTx[keepID,]
  return(myfusionTx)
}



convertChrPos <- function(txname,txpos,exonInfo=NULL,readLen=100,txExonMat=NULL){
  ##### Converting trancript positions to chromosome positions. Be aware of the direction of strands and read mapping
  if (is.null(txExonMat)) txExonMat=exonInfo[which(exonInfo$TXNAME==txname),]
  
  if (txExonMat$TXSTRAND[1]=="+") {
    txExonMat=txExonMat[order(txExonMat$EXONSTART),]#sort exons by increasing order for forward strand
    txlen=sum(txExonMat$EXONEND-txExonMat$EXONSTART+1);
    exonlen=txExonMat$EXONEND-txExonMat$EXONSTART+1
    exonlenCumSum=cumsum(exonlen)
    chrPos=NULL
    for (i in 1:length(txpos)){
      exID=which(exonlenCumSum>txpos[i])[1]
      if (!is.na(exID)) posOffset=exonlen[exID]-(exonlenCumSum[exID]-txpos[i]) 
      else{
        #if position outside the range, use the last exon
        exID=length(exonlen)
        posOffset=exonlen[exID]+(txpos[i]-exonlenCumSum[exID]) 
      }
      chrPos=c(chrPos,txExonMat$EXONSTART[exID]+posOffset)
    }
    return(chrPos)
  }else{ 
    txExonMat=txExonMat[order(txExonMat$EXONEND, decreasing=TRUE),]#sort exons by decreasing order for reverse strand
    txlen=sum(txExonMat$EXONEND-txExonMat$EXONSTART+1);
    exonlen=txExonMat$EXONEND-txExonMat$EXONSTART+1
    exonlenCumSum=cumsum(exonlen)
    chrPos=NULL
    for (i in 1:length(txpos)){
      exID=which(exonlenCumSum>txpos[i])[1]
      if (!is.na(exID)) posOffset=exonlenCumSum[exID]-txpos[i]-1
      else{
        exID=length(exonlen)
        posOffset=exonlenCumSum[exID]-txpos[i]
      }
      chrPos=c(chrPos,txExonMat$EXONSTART[exID]+posOffset)
    }
    return(chrPos)
  }
}


convertGenePos <- function(txname,txpos,genename,exonInfo,readLen=100){
  ##### Converting trancript positions to gene positions. Be aware of the direction of strands and read mapping
  mygeneEx=exonInfo[exonInfo$GENEID==genename,]
  
  if (mygeneEx$TXSTRAND[1]=="+") {
    mygeneEx=mygeneEx[order(mygeneEx$EXONID),] #sort exons by increasing order for forward strand
    brStart=unique(mygeneEx$EXONSTART)
    brEnd=unique(mygeneEx$EXONEND)
    exStartStatus=sapply(brStart, function(x) sum((x-mygeneEx$EXONSTART)/1e8*(x-mygeneEx$EXONEND)<0))
    exEndStatus=sapply(brEnd, function(x) sum((x-mygeneEx$EXONSTART)/1e8*(x-mygeneEx$EXONEND)<0))
    brStart=brStart[exStartStatus==0] #none-overlap start
    brEnd=brEnd[exEndStatus==0]#none-overlap end
    brCumsum=cumsum(brEnd-brStart+1)
    mytx=mygeneEx[mygeneEx$TXNAME==txname,]
    mytxChrPos=min(mytx$EXONSTART)
    brID=which((mytxChrPos-brStart)/1e8*(mytxChrPos-brEnd)<=0) #which segment the transcript belongs to
    mytxGenePos=ifelse(brID==1,0,brCumsum[brID-1])+mytxChrPos-brStart[brID] #position of the tx at the gene
    genePos=min(mytxGenePos)+txpos 
    return(genePos)
  }else{
    mygeneEx=mygeneEx[order(mygeneEx$EXONID, decreasing=TRUE),] #sort exons by decreasing order for reverse strand
    brStart=unique(mygeneEx$EXONSTART)
    brEnd=unique(mygeneEx$EXONEND)
    exStartStatus=sapply(brStart, function(x) sum((x-mygeneEx$EXONSTART)/1e8*(x-mygeneEx$EXONEND)<0))
    exEndStatus=sapply(brEnd, function(x) sum((x-mygeneEx$EXONSTART)/1e8*(x-mygeneEx$EXONEND)<0))
    brStart=brStart[exStartStatus==0] #none-overlap start
    brEnd=brEnd[exEndStatus==0] #none-overlap end
    brCumsum=cumsum(brEnd-brStart+1)
    mytx=mygeneEx[mygeneEx$TXNAME==txname,]
    mytxlen=sum(mytx$EXONEND-mytx$EXONSTART+1)
    mytxChrPos=min(mytx$EXONSTART)
    brID=which((mytxChrPos-brStart)/1e8*(mytxChrPos-brEnd)<=0)#which segment the transcript belongs to
    mytxGenePos=max(brCumsum)-(brCumsum[brID]-(mytxChrPos-brStart[brID])) #position of the tx at the gene
    genePos=min(mytxGenePos) + mytxlen - txpos #This should be the leftmost position of the gene
    return(genePos)
  }
}


convertChrPosGenePos<- function (chrPos,genename=NULL,exonInfo=NULL,geneExonMat=NULL){
  if (is.null(geneExonMat)) geneExonMat=exonInfo[exonInfo$GENEID==genename,]
  if (geneExonMat$TXSTRAND[1]=="+") {
    geneExonMat=geneExonMat[order(geneExonMat$EXONID),] #sort exons by increasing order for forward strand
    brStart=unique(geneExonMat$EXONSTART)
    brEnd=unique(geneExonMat$EXONEND)
    exStartStatus=sapply(brStart, function(x) sum((x-geneExonMat$EXONSTART)/1e8*(x-geneExonMat$EXONEND)<0))
    exEndStatus=sapply(brEnd, function(x) sum((x-geneExonMat$EXONSTART)/1e8*(x-geneExonMat$EXONEND)<0))
    brStart=brStart[exStartStatus==0] #none-overlap start
    brEnd=brEnd[exEndStatus==0]#none-overlap end
    brCumsum=cumsum(brEnd-brStart+1)
    genePos=NULL;
    for (i in 1:length(chrPos)){
      brID=which((chrPos[i]-brStart)/1e8*(chrPos[i]-brEnd)<=0) #which segment the chrPos belongs to
      myGenePos=ifelse(brID==1,0,brCumsum[brID-1])+chrPos[i]-brStart[brID] #position of the chrPos at the gene
      genePos=c(genePos,min(myGenePos))
    }
    return(genePos)
  }else{
    geneExonMat=geneExonMat[order(geneExonMat$EXONID, decreasing=TRUE),] #sort exons by decreasing order for reverse strand
    brStart=unique(geneExonMat$EXONSTART)
    brEnd=unique(geneExonMat$EXONEND)
    exStartStatus=sapply(brStart, function(x) sum((x-geneExonMat$EXONSTART)/1e8*(x-geneExonMat$EXONEND)<0))
    exEndStatus=sapply(brEnd, function(x) sum((x-geneExonMat$EXONSTART)/1e8*(x-geneExonMat$EXONEND)<0))
    brStart=brStart[exStartStatus==0] #none-overlap start
    brEnd=brEnd[exEndStatus==0] #none-overlap end
    brCumsum=cumsum(brEnd-brStart+1)
    genePos=NULL;
    for (i in 1:length(chrPos)){
      brID=which((chrPos[i]-brStart)/1e8*(chrPos[i]-brEnd)<=0)#which segment the chrPos belongs to
      myGenePos=max(brCumsum)-(brCumsum[brID]-(chrPos[i]-brStart[brID])) #position of the chrPos at the gene
      genePos=c(genePos,min(myGenePos)) #This should be the leftmost position of the gene
    }
    return(genePos)
  }
  
}


detectExonID <- function(chrPos,txname,exonInfo){
  myTxEx=exonInfo[exonInfo$TXNAME==txname,]
  exID=lapply(chrPos,function(x){ 
    myExid=which((myTxEx$EXONSTART-x)/1e8*(myTxEx$EXONEND-x)<=0)
    return(myTxEx$EXONID[myExid])
  })
  exID=unique(unlist(exID))
  return(exID)
}

detectExonID2 <- function(chrPos,seqLen,txname,exonInfo){
  myTxEx=exonInfo[exonInfo$TXNAME==txname,]
  if (myTxEx$TXSTRAND[1]=="+") chrPos2=chrPos+seqLen else chrPos2=chrPos-seqLen
  
  exID=lapply(cbind(chrPos,chrPos2),function(x){ 
    myExid1=which((myTxEx$EXONSTART-x[1])/1e8*(myTxEx$EXONEND-x[1])<=0)
    myExid2=which((myTxEx$EXONSTART-x[2])/1e8*(myTxEx$EXONEND-x[2])<=0)
    myExid=unique(c(myExid1,myExid2))
    return(myTxEx$EXONID[myExid])
  })
  exID=unique(unlist(exID))
  return(exID)
}


detectExonID3 <- function(txname,txpos,seqLen, exonInfo){
  txExonMat=exonInfo[which(exonInfo$TXNAME==txname),]
  
  if (txExonMat$TXSTRAND[1]=="+") {
    txExonMat=txExonMat[order(txExonMat$EXONSTART),]#sort exons by increasing order for forward strand
    txlen=sum(txExonMat$EXONEND-txExonMat$EXONSTART+1);
    exonlen=txExonMat$EXONEND-txExonMat$EXONSTART+1
    exonlenCumSum=cumsum(exonlen)
    exList=NULL
    for (i in 1:length(txpos)){
      exID=which(exonlenCumSum>txpos[i])[1]
      if (!is.na(exID)){
        exID2=which(exonlenCumSum>txpos[i]+seqLen[i])[1]
        if (is.na(exID2)) exID2=length(exonlen)
        exList=c(exList,exID:exID2)
        
      }
      else{
        #if position outside the range, use the last exon
        exID=length(exonlen)
        exList=c(exList,exID)
      }
    }
    return(txExonMat$EXONID[unique(exList)])
  }else{ 
    txExonMat=txExonMat[order(txExonMat$EXONEND, decreasing=TRUE),]#sort exons by decreasing order for reverse strand
    txlen=sum(txExonMat$EXONEND-txExonMat$EXONSTART+1);
    exonlen=txExonMat$EXONEND-txExonMat$EXONSTART+1
    exonlenCumSum=cumsum(exonlen)
    exList=NULL
    for (i in 1:length(txpos)){
      #txpos[i]=txlen-txpos[i]-readLen #convert pos3 to pos5
      exID=which(exonlenCumSum>txpos[i])[1]
      if (!is.na(exID)){
        exID2=which(exonlenCumSum>txpos[i]+seqLen[i])[1]
        if (is.na(exID2)) exID2=length(exonlen)
        exList=c(exList,exID:exID2)
      }
      else{
        #if position outside the range, use the last exon
        exID=length(exonlen)
        exList=c(exList,exID)
      }
    }
    return(txExonMat$EXONID[unique(exList)])
  }
}



getGeneLen <-function(genename=NULL,exonInfo=NULL,geneExonMat=NULL){
  if (is.null(geneExonMat)) geneExonMat=exonInfo[exonInfo$GENEID==genename,]
  if (geneExonMat$TXSTRAND[1]=="+") {
    geneExonMat=geneExonMat[order(geneExonMat$EXONID),] #sort exons by increasing order for forward strand
    brStart=unique(geneExonMat$EXONSTART)
    brEnd=unique(geneExonMat$EXONEND)
    exStartStatus=sapply(brStart, function(x) sum((x-geneExonMat$EXONSTART)/1e8*(x-geneExonMat$EXONEND)<0))
    exEndStatus=sapply(brEnd, function(x) sum((x-geneExonMat$EXONSTART)/1e8*(x-geneExonMat$EXONEND)<0))
    brStart=brStart[exStartStatus==0] #none-overlap start
    brEnd=brEnd[exEndStatus==0]#none-overlap end
    brCumsum=cumsum(brEnd-brStart+1)
    return(max(brCumsum))
  }else{
    geneExonMat=geneExonMat[order(geneExonMat$EXONID, decreasing=TRUE),] #sort exons by decreasing order for reverse strand
    brStart=unique(geneExonMat$EXONSTART)
    brEnd=unique(geneExonMat$EXONEND)
    exStartStatus=sapply(brStart, function(x) sum((x-geneExonMat$EXONSTART)/1e8*(x-geneExonMat$EXONEND)<0))
    exEndStatus=sapply(brEnd, function(x) sum((x-geneExonMat$EXONSTART)/1e8*(x-geneExonMat$EXONEND)<0))
    brStart=brStart[exStartStatus==0]
    brEnd=brEnd[exEndStatus==0]
    brCumsum=cumsum(brEnd-brStart+1)
    return(max(brCumsum))
  }
}



# Filtering by dupGene
filterBydupGene <- function(myfusionTx,dupGene.thres=4){
  dupGene1=table(myfusionTx$gene1)
  dupGene2=table(myfusionTx$gene2)
  myfusionTx=cbind(myfusionTx,as.integer(dupGene1[match(as.character(myfusionTx$gene1),names(dupGene1))]))
  myfusionTx=cbind(myfusionTx,as.integer(dupGene2[match(as.character(myfusionTx$gene2),names(dupGene2))]))
  colnames(myfusionTx)[c(ncol(myfusionTx)-1,ncol(myfusionTx))]=c("dupGene1","dupGene2")
  keepID=which(myfusionTx$dupGene1 <= dupGene.thres & myfusionTx$dupGene2<=dupGene.thres)
  myfusionTx=myfusionTx[keepID,]
  return(myfusionTx)
}

computeDupGene <- function(myfusionTx,dupGene.thres=4){
  dupGene1=table(myfusionTx$gene1)
  dupGene2=table(myfusionTx$gene2)
  myfusionTx=cbind(myfusionTx,as.integer(dupGene1[match(as.character(myfusionTx$gene1),names(dupGene1))]))
  myfusionTx=cbind(myfusionTx,as.integer(dupGene2[match(as.character(myfusionTx$gene2),names(dupGene2))]))
  colnames(myfusionTx)[c(ncol(myfusionTx)-1,ncol(myfusionTx))]=c("dupGene1","dupGene2")
  
  if (dupGene.thres<0) return(myfusionTx[,c(ncol(myfusionTx)-1,ncol(myfusionTx))])
  
  keepID=which(myfusionTx$dupGene1 <= dupGene.thres & myfusionTx$dupGene2<=dupGene.thres)
  myfusionTx=myfusionTx[keepID,]
  return(myfusionTx)
}


computeReversedFusionCount <- function(myfusionTx){
  #similar to filterByReverseStrand, however, we consider only gene-level with real count
  matchID=match(myfusionTx$name12, myfusionTx$name21)
  name21Count=myfusionTx$supportCount[matchID]
  name21Count[is.na(name21Count)]=0
  return(name21Count)
}


computeGeneDistance <- function(myfusionTx,anntxdb,minGeneDist=1e5){
  res1=select(anntxdb, keys=as.character(myfusionTx$gene1), columns=c("TXSTART","TXEND"), keytype = "GENEID")
  res2=select(anntxdb, keys=as.character(myfusionTx$gene2), columns=c("TXSTART","TXEND"), keytype = "GENEID")
  
  minRes1=tapply(res1$TXSTART,as.character(res1$GENEID),min)
  maxRes1=tapply(res1$TXEND,as.character(res1$GENEID),max)
  minRes1=minRes1[match(as.character(myfusionTx$gene1),names(minRes1))]
  maxRes1=maxRes1[match(as.character(myfusionTx$gene1),names(maxRes1))]
  
  minRes2=tapply(res2$TXSTART,as.character(res2$GENEID),min)
  maxRes2=tapply(res2$TXEND,as.character(res2$GENEID),max)
  minRes2=minRes2[match(as.character(myfusionTx$gene2),names(minRes2))]
  maxRes2=maxRes2[match(as.character(myfusionTx$gene2),names(maxRes2))]
  
  geneDist=unlist(apply(cbind(abs(minRes1-maxRes2),abs(minRes1-minRes2),abs(maxRes1-maxRes2),abs(maxRes1-minRes2)),1,min))
  if (minGeneDist < 0) return(geneDist);
  
  myfusionTx$geneDist=geneDist
  myfusionTx=myfusionTx[myfusionTx$geneDist>=minGeneDist,]
  return(myfusionTx)
}



estimateCountEM <- function(alp0,alpOut,feq,fgeFeqMap,itNum=2,alpDiff.thres=0.01){
  ##### Estimating count from fusion-equivalence classes using EM
  isConvergered=FALSE
  alpIn=alp0
  for (k in 1:itNum){
    for (i in 1:length(feq)){
      feqID=i;
      feqCount=feq[feqID]
      geneID=unlist(fgeFeqMap[feqID])
      if (length(geneID) > 1){
        weights=rep(1/length(geneID),length(geneID))
        v=alpIn[geneID] * weights;
        denom=sum(v)
        invt_denom = feqCount / denom;
        alpOut[geneID]=alpOut[geneID]+v*invt_denom
      }else{
        alpOut[geneID]=alpOut[geneID]+feqCount
      }
    }
    alpDiff = abs(alpIn - alpOut) / alpOut;
    if (sum(alpDiff > alpDiff.thres) == 0){
      isConvergered=TRUE
      k=itNum+1
    }
    alpIn=alpOut
    alpOut=rep(0,length(alpOut))
  }
  res=list(correctedCount=alpIn,isConvergered=isConvergered)
  return(res)
}


convertReverseComplement<-function(DNAseq){
  DNAarr=unlist(strsplit(DNAseq,""))
  #reverse
  DNAarr=rev(DNAarr)
  #complement
  Aid=which(DNAarr=="A")
  Tid=which(DNAarr=="T")
  Gid=which(DNAarr=="G")
  Cid=which(DNAarr=="C")
  DNAarr[Aid]="T"
  DNAarr[Tid]="A"
  DNAarr[Gid]="C"
  DNAarr[Cid]="G"
  #result
  DNAseqRc=paste(DNAarr,collapse = "")
 return(DNAseqRc) 
}



matchFusionReads <- function(inPath, readStrands, fastaOut, junctInfo, splitReads, fusionCandidate,fsizeLadder){
  ##### do "pseudo de novo" alignment based on the mapped positions of reads and match reads between 5 prime and 3 prime sites
  fusionName=fusionCandidate$fusionName
  fastaDat2=fastaDat1=list();
  fastaSplit2=fastaSplit1=list();
  frfiles=list.files(inPath,paste(readStrands,"_fusionMappedReadsChunk_*",sep=""))
  for (i in 1:length(frfiles)){
    #read fasta files of mapped reads
    ftag=rev(strsplit(strsplit(frfiles[i],"\\.")[[1]][1],"_")[[1]])[1]
    con <- file(paste(inPath,"/",readStrands,"_fastaseq_",ftag,"_1.fa",sep=""), "r", blocking = FALSE)
    mydata=readLines(con)
    close(con)
    fastaDat1[[frfiles[i]]]=mydata
    
    con <- file(paste(inPath,"/",readStrands,"_fastaseq_",ftag,"_2.fa",sep=""), "r", blocking = FALSE)
    mydata=readLines(con)
    close(con)
    fastaDat2[[frfiles[i]]]=mydata
    
    #read fasta files of split reads
    con <- file(paste(inPath,"/","splitRead_",ftag,"_1.fa",sep=""), "r", blocking = FALSE)
    mydata=readLines(con)
    close(con)
    fastaSplit1[[frfiles[i]]]=mydata
    
    con <- file(paste(inPath,"/","splitRead_",ftag,"_2.fa",sep=""), "r", blocking = FALSE)
    mydata=readLines(con)
    close(con)
    fastaSplit2[[frfiles[i]]]=mydata
    
  }
  
  fastaSplit1=unlist(fastaSplit1)
  fastaSplit2=unlist(fastaSplit2)
  
  
  ###### de novo assembly by alignment
  matchInfo=list();
  consensusEntropy=consensusErr=consensusProp=NULL
  for (i in 1:nrow(fusionCandidate)){
    fName=fusionCandidate$fusionName[i]
    res=junctInfo[[fName]]
    myreadID=res$readID
    rID=unlist(myreadID)
    fID=unlist(lapply(rID, function(x) which(fsizeLadder>=x)[1]))
    rID.adj=unlist(lapply(seq_along(rID), function(x) if (fID[x] > 1) return(rID[x]-fsizeLadder[fID[x]-1]) else return(rID[x])))
    
    if (readStrands=="FR" ||readStrands=="FF"){ #read1 for 5'  and read 2 for 3'
      if (fusionCandidate$strand1[i]=="+"){
        offset1=-unlist(res$readL.seqLen)+1-unlist(res$readL.Pos)+unlist(res$readL.ReadPos)
        mypos1=unlist(res$adj.readL.GenePos)+offset1
      }else{
        offset1=unlist(res$readL.seqLen)-1+unlist(res$readL.Pos)-unlist(res$readL.ReadPos)
        mypos1=unlist(res$adj.readL.GenePos)+offset1
      }
      if (fusionCandidate$strand2[i]=="+"){
        offset2=-unlist(res$readR.Pos)+unlist(res$readR.ReadPos)
        mypos2=unlist(res$adj.readR.GenePos)+offset2
      } else{
        offset2=unlist(res$readR.Pos)-unlist(res$readR.ReadPos)
        mypos2=unlist(res$adj.readR.GenePos)+offset2
      }
      if (fusionCandidate$strand1[i]=="-"){
        mypos1=-mypos1
        offset1=-offset1
      }
      if (fusionCandidate$strand2[i]=="-"){
        mypos2=-mypos2
        offset2=-offset2
      }
      
      mypos1=mypos1-min(mypos1)
      mypos2=mypos2-min(mypos2)
      
      mySeq1=mySeq2=NULL
      mySeq1N=mySeq2N=NULL
      for (j in 1:length(rID.adj)){
        mySeq1=c(mySeq1,fastaDat1[[fID[j]]][rID.adj[j]*2])
        mySeq2=c(mySeq2,fastaDat2[[fID[j]]][rID.adj[j]*2])
        
        mySeq1N=c(mySeq1N,paste(paste(rep("N",mypos1[j]),collapse = ""),mySeq1[j],sep=""))
        
        mySeqRc=convertReverseComplement(mySeq2[j]) #make sure read2 is reverse complement
        mySeq2N=c(mySeq2N,paste(paste(rep("N",mypos2[j]),collapse = ""),mySeqRc,sep=""))
      }
      #Sometimes, the alignment is not optimal, we might need a re-alignment step here
      estBr5=max(mypos1-offset1)
      estBr3=min(mypos2-offset2)
      mySeq5N=mySeq1N
      mySeq3N=mySeq2N
      #####get split reads
      myreadID=which(splitReads$name12==fName)
      myheader=as.character(splitReads$header[myreadID])
      rID.adj=match(myheader,fastaSplit1)
      mySplit5N=fastaSplit1[rID.adj+1]
      mySplit3N=fastaSplit2[rID.adj+1]
      
    } else{ #if readStrands=RF, UN or RR, so read1 for 3'  and read 2 for 5'
      
      if (fusionCandidate$strand1[i]=="+"){
        #offset1=-unlist(res$readL.seqLen)+1-unlist(res$readL.Pos)+unlist(res$readL.ReadPos)
        offset1=-unlist(res$readL.Pos)+unlist(res$readL.ReadPos)
        mypos1=unlist(res$adj.readL.GenePos)+offset1
      }else{
        #offset1=unlist(res$readL.seqLen)-1+unlist(res$readL.Pos)-unlist(res$readL.ReadPos)
        offset1=unlist(res$readL.Pos)-unlist(res$readL.ReadPos)
        mypos1=unlist(res$adj.readL.GenePos)+offset1
      }
      
      if (fusionCandidate$strand2[i]=="+"){
        offset2=-unlist(res$readR.seqLen)+1-unlist(res$readR.Pos)+unlist(res$readR.ReadPos)
        mypos2=unlist(res$adj.readR.GenePos)+offset2
      }else{
        offset2=unlist(res$readR.seqLen)-1+unlist(res$readR.Pos)-unlist(res$readR.ReadPos)
        mypos2=unlist(res$adj.readR.GenePos)+offset2
      }
      
      if (fusionCandidate$strand1[i]=="-"){
        mypos1=-mypos1
        offset1=-offset1
      }
      if (fusionCandidate$strand2[i]=="-"){
        mypos2=-mypos2
        offset2=-offset2
      }
      
      mypos1=mypos1-min(mypos1)
      mypos2=mypos2-min(mypos2)
      
      mySeq1=mySeq2=NULL
      mySeq1N=mySeq2N=NULL
      for (j in 1:length(rID.adj)){
        mySeq1=c(mySeq1,fastaDat1[[fID[j]]][rID.adj[j]*2])
        mySeq2=c(mySeq2,fastaDat2[[fID[j]]][rID.adj[j]*2])
        mySeqRc=convertReverseComplement(mySeq1[j]) #make sure read1 is reverse complement
        mySeq1N=c(mySeq1N,paste(paste(rep("N",mypos1[j]),collapse = ""),mySeqRc,sep=""))
        
        mySeq2N=c(mySeq2N,paste(paste(rep("N",mypos2[j]),collapse = ""),mySeq2[j],sep=""))
        
      }
      #Sometimes, the alignment is not optimal, we might need a re-alignment step here
      estBr5=max(mypos2-offset2)
      estBr3=min(mypos1-offset1)
      mySeq5N=mySeq2N
      mySeq3N=mySeq1N
      
      #####get split reads
      myreadID=which(splitReads$name12==fName)
      myheader=as.character(splitReads$header[myreadID])
      rID.adj=match(myheader,fastaSplit1)
      mySplit5N=fastaSplit2[rID.adj+1]
      mySplit3N=fastaSplit1[rID.adj+1]
      
    }
    
    
    ##### do alignment
    
    myseq5len=max(nchar(mySeq5N))
    myseq3len=max(nchar(mySeq3N))
    res5=sapply(mySeq5N,function(x) strsplit(x,""))
    res3=sapply(mySeq3N,function(x) strsplit(x,""))
    
    res5=lapply(res5,function(x) x=c(x,rep("N",myseq5len-length(x))))
    res3=lapply(res3,function(x) x=c(x,rep("N",myseq3len-length(x))))
    #get nucleotide frequencies
    nuc <- c("A","T","G","C","N")
    res=do.call(rbind,res5)
    nucFreq5=NULL
    for (j in 1:ncol(res)){
      nucFreq5=rbind(nucFreq5,table(c(nuc,res[,j])))
    }
    nucFreq5=nucFreq5-1
    
    res=do.call(rbind,res3)
    nucFreq3=NULL
    for (j in 1:ncol(res)){
      nucFreq3=rbind(nucFreq3,table(c(nuc,res[,j])))
    }
    nucFreq3=nucFreq3-1
    
    nucFreq5=nucFreq5[,-4]
    nucFreq3=nucFreq3[,-4]
    #find the break points
    # start5=estBr5-1
    # start3=estBr3+1
    start5=myseq5len
    start3=myseq3len
    
    fcErr=fcProp=NULL
    entropy3=entropy5=entropyInter=NULL
    for (j in 10:start3){ #at least 10 bases (1/1e6 by chance) overlapped between two transcripts
      newstart3=j
      #get new match
      inter5StarPos=start5-(newstart3-1)
      if (inter5StarPos < 1) inter5StarPos=1
      inter5=nucFreq5[inter5StarPos:myseq5len,]
      
      inter3EndPos=newstart3+myseq5len-start5 # we fix start5=myseq5len, so it will add zero
      if (inter3EndPos > myseq3len) inter3EndPos=myseq3len
      inter3=nucFreq3[1:inter3EndPos,]
      
      if (nrow(inter5)!=nrow(inter3)) break();
      
      consProp5=apply(inter5,1,max)/rowSums(inter5)
      consProp3=apply(inter3,1,max)/rowSums(inter3)
      consProp5[is.na(consProp5)]=1
      consProp3[is.na(consProp3)]=1
      consPropInter=apply(inter3+inter5,1,max)/rowSums(inter3+inter5)
      consPropInter[is.na(consPropInter)]=1e-9
      
      fcProp=c(fcProp,(sum(consProp5>consPropInter)+sum(consProp3>consPropInter))/(2*length(consPropInter)))
      fcErr=c(fcErr,(mean(consProp5-consPropInter)+mean(consProp3-consPropInter))/2)
      #      fcProp=c(fcProp,2*mean(consPropInter)/(mean(consProp5)+mean(consProp3)))
      entropyInter=c(entropyInter,-sum(consPropInter*log(consPropInter)))
      entropy3=c(entropy3,-sum(consProp3*log(consProp3)))
      entropy5=c(entropy5,-sum(consProp5*log(consProp5)))
    }
    
 
    res=list()
    res$estBr5=estBr5
    res$estBr3=estBr3
    res$myseq5len=myseq5len
    res$myseq3len=myseq3len
    res$mySeq5N=mySeq5N
    res$mySeq3N=mySeq3N
    res$mySplit5N=mySplit5N
    res$mySplit3N=mySplit3N
    res$nucFreq5=nucFreq5
    res$nucFreq3=nucFreq3
    res$fcProp=fcProp
    res$fcErr=fcErr
    res$entropyInter=entropyInter
    res$entropy3=entropy3
    res$entropy5=entropy5
    
    matchInfo[[fName]]=res
    
    
    consensusProp=c(consensusProp,min(fcProp))
    consensusErr=c(consensusErr,min(fcErr))
    consensusEntropy=c(consensusEntropy,min(entropyInter))
    
  }
  
  fusionCandidate$consensusEntropy=consensusEntropy
  fusionCandidate$consensusProp=consensusProp
  fusionCandidate$consensusErr=consensusErr
  fusionCandidate$consensusScore=fusionCandidate$consensusErr*fusionCandidate$consensusProp

  return(list(fusionCandidate=fusionCandidate,matchInfo=matchInfo))
}



exportMappedFusionReads <- function(inPath, readStrands, fastaOut, junctInfo, fusionName,fsizeLadder){
  ##### export mapped reads of fusion genes to files
  fastaDat2=fastaDat1=list();
  frfiles=list.files(inPath,paste(readStrands,"_fusionMappedReadsChunk_*",sep=""))
  for (i in 1:length(frfiles)){
    #read fasta files
    ftag=rev(strsplit(strsplit(frfiles[i],"\\.")[[1]][1],"_")[[1]])[1]
    con <- file(paste(inPath,"/",readStrands,"_fastaseq_",ftag,"_1.fa",sep=""), "r", blocking = FALSE)
    mydata=readLines(con)
    close(con)
    fastaDat1[[frfiles[i]]]=mydata
    
    con <- file(paste(inPath,"/",readStrands,"_fastaseq_",ftag,"_2.fa",sep=""), "r", blocking = FALSE)
    mydata=readLines(con)
    close(con)
    fastaDat2[[frfiles[i]]]=mydata
  }
  
  if (is.null(fastaOut) || is.na(fastaOut)) fastaOut=""
  con1 <- file(paste(fastaOut,readStrands,"_fusionReads_1.fa",sep=""), "w", blocking = FALSE)
  con2 <- file(paste(fastaOut,readStrands,"_fusionReads_2.fa",sep=""), "w", blocking = FALSE)
  
  ###########  
  #detect read positions in fasta files
  for (fName in fusionName){
    res=junctInfo[[fName]]
    if (length(res)>0){
      myreadID=res$readID
      rID=unlist(myreadID)
      rID=sort(rID)
      fID=unlist(lapply(rID, function(x) which(fsizeLadder>=x)[1]))
      rID.adj=unlist(lapply(seq_along(rID), function(x) if (fID[x] > 1) return(rID[x]-fsizeLadder[fID[x]-1]) else return(rID[x])))
      #export reads to file
      for (j in 1:length(rID.adj)){
        readHead=paste(">FuSeq__MR__",unlist(res$readL.ReadPos)[j],"__",unlist(res$readR.ReadPos)[j],"__",fName," /1",sep="")
        writeLines(readHead,con1)
        writeLines(fastaDat1[[fID[j]]][rID.adj[j]*2],con1)
        
        readHead=paste(">FuSeq__MR__",unlist(res$readL.ReadPos)[j],"__",unlist(res$readR.ReadPos)[j],"__",fName," /2",sep="")
        writeLines(readHead,con2)
        writeLines(fastaDat2[[fID[j]]][rID.adj[j]*2],con2)
      }
    }
  }
  close(con1)
  close(con2)
}



exportSplitFusionReads <- function(inPath, readStrands, fastaOut, splitReads, fusionName){
  ##### export read paids of the split reads of fusion genes to files
  fastaDat2=fastaDat1=list();
  frfiles=list.files(inPath,paste("splitReadInfo_*",sep=""))
  for (i in 1:length(frfiles)){
    #read fasta files
    ftag=rev(strsplit(strsplit(frfiles[i],"\\.")[[1]][1],"_")[[1]])[1]
    con <- file(paste(inPath,"/","splitRead_",ftag,"_1.fa",sep=""), "r", blocking = FALSE)
    mydata=readLines(con)
    close(con)
    fastaDat1[[frfiles[i]]]=mydata
    
    con <- file(paste(inPath,"/","splitRead_",ftag,"_2.fa",sep=""), "r", blocking = FALSE)
    mydata=readLines(con)
    close(con)
    fastaDat2[[frfiles[i]]]=mydata
  }
  
  fastaDat1=unlist(fastaDat1)
  fastaDat2=unlist(fastaDat2)
  
  if (is.null(fastaOut) || is.na(fastaOut)) fastaOut=""
  con1 <- file(paste(fastaOut,readStrands,"_fusionReads_1.fa",sep=""), "a", blocking = FALSE)
  con2 <- file(paste(fastaOut,readStrands,"_fusionReads_2.fa",sep=""), "a", blocking = FALSE)
  
  
  #detect read positions in fasta files
  for (fName in fusionName){
    myreadID=which(splitReads$name12==fName)
    if (length(myreadID)>0){
      myheader=as.character(splitReads$header[myreadID])
      rID.adj=match(myheader,fastaDat1)
      #export reads to file
      for (j in rID.adj){
        readHead=paste(">FuSeq__SR__",splitReads$front_hitpos[j],"__",splitReads$back_hitpos[j],"__",fName," /1",sep="")
        writeLines(readHead,con1)
        writeLines(fastaDat1[j+1],con1)
        
        readHead=paste(">FuSeq__SR__",splitReads$front_hitpos[j],"__",splitReads$back_hitpos[j],"__",fName," /2",sep="")
        writeLines(readHead,con2)
        writeLines(fastaDat2[j+1],con2)
      }
    }
  }
  close(con1)
  close(con2)
}


testFtxlen <- function(mu, sig, r, ftxlen, kmerlen,fragDist,M=10000){
  # to get a distribution of one side fragments with the limitation of k-mer length 
  x = runif(M, min=0, max= min(ftxlen,mu+2*sig)) 
  len = sample(fragDist[,1], M, replace=TRUE, prob = fragDist[,2])
  #limit x by the kmerlen
  keepID=which(x>=kmerlen)
  x=x[keepID]
  len=len[keepID]
  return(list(x=x, len=len))
}



refineJunctionBreak <-function(junctBr, anntxdb, fragmentInfo, readStrands, fragDist){
  ##########################
  #### revisit junction breaks to remove outliers and detect splicing sites
  ##########################
  exonInfo=select(anntxdb, keys=unique(c(as.character(junctBr$myFusionFinal$gene1),as.character(junctBr$myFusionFinal$gene2))), columns=c("TXNAME","EXONID","EXONSTART","EXONEND","TXSTRAND"), keytype = "GENEID")
  
  junctBr$myFusionFinal$adjFragmentNum=junctBr$myFusionFinal$outlierNum=NA
  junctBr$myFusionFinal$ftxLenProp=junctBr$myFusionFinal$ftxMedianLen=junctBr$myFusionFinal$ftxMeanLen=NA
  junctBr$myFusionFinal$ftx5LenSd=junctBr$myFusionFinal$ftx5LenMean=junctBr$myFusionFinal$ftx3LenSd=junctBr$myFusionFinal$ftx3LenMean=NA
  
  junctBr$myFusionFinal$crt.tx5GapMean=junctBr$myFusionFinal$crt.tx3GapMean=junctBr$myFusionFinal$crt.tx5GapSd=junctBr$myFusionFinal$crt.tx3GapSd=junctBr$myFusionFinal$crt.tx5LenMean=junctBr$myFusionFinal$crt.tx3LenMean=junctBr$myFusionFinal$crt.tx5LenSd=junctBr$myFusionFinal$crt.tx3LenSd=NA  
  for (i in 1:nrow(junctBr$myFusionFinal)){
    #which(junctBr$myFusionFinal$fusionName=="ENSG00000076864-ENSG00000138079")
    fgeName=as.character(junctBr$myFusionFinal$name12)[i]
    gene1=as.character(junctBr$myFusionFinal$gene1[i])
    gene2=as.character(junctBr$myFusionFinal$gene2[i])
    
    res=junctBr$junctInfo[[i]]
    
    #detect outliers
    adjPosL=unlist(res$adj.readL.GenePos)
    adjPosR=unlist(res$adj.readR.GenePos)
    rmID=NULL
    
    #need empirical fragment length distribution here
    flenRef=rep(fragDist[,1],fragDist[,2])
    flenRef=flenRef-mean(flenRef)
    
    ###### use qqplot to detect outliers: build a line by y=bx+a and compute the distance of value to the line
    # flenRef=flenRef[flenRef<quantile(flenRef,0.95) & flenRef>quantile(flenRef,0.05)]
    # quanValRef=quantile(flenRef,probs=c(0.25,0.75))
    # 
    # #left
    # adjPosL=adjPosL-median(adjPosL)
    # myorder=order(adjPosL,decreasing = FALSE)
    # adjPosL=adjPosL[myorder]
    # quanVal=quantile(adjPosL,probs=c(0.25,0.75))
    # qqpos = qqplot(flenRef,adjPosL,plot.it = FALSE)
    # bval=diff(quanVal)/diff(quanValRef)
    # a=quanVal[1]-bval*quanValRef[1]
    # ref.val = a+ bval*qqpos$x
    # #points(qqpos$x, ref.val, col='red')
    # res.val = qqpos$y - ref.val
    # sig = diff(quanVal)/(qnorm(0.75)-qnorm(0.25))
    # stdres = res.val/sig
    # outlierID=which(stdres< -2)
    # #adjPosL[outlierID]
    # outlierID=myorder[outlierID]
    # #get back the original values of adjPosL
    # adjPosL=unlist(res$adj.readL.GenePos)
    # #adjPosL=adjPosL[myorder[order(myorder)]]
    # #adjPosL[outlierID]
    # rmID=c(rmID,outlierID)
    # 
    # #right
    # adjPosR=adjPosR-median(adjPosR)
    # myorder=order(adjPosR,decreasing = FALSE)
    # adjPosR=adjPosR[myorder]
    # quanVal=quantile(adjPosR,probs=c(0.25,0.75))
    # qqpos = qqplot(flenRef,adjPosR,plot.it = FALSE)
    # bval=diff(quanVal)/diff(quanValRef)
    # a=quanVal[1]-bval*quanValRef[1]
    # ref.val = a+ bval*qqpos$x
    # #points(qqpos$x, ref.val, col='red')
    # res.val = qqpos$y - ref.val
    # sig = diff(quanVal)/(qnorm(0.75)-qnorm(0.25))
    # stdres = res.val/sig
    # outlierID=which(stdres< -2)
    # #adjPosR[outlierID]
    # outlierID=myorder[outlierID]
    # #get back the original values of adjPosL
    # adjPosR=unlist(res$adj.readR.GenePos)
    # #adjPosR=adjPosR[myorder[order(myorder)]]
    # #adjPosR[outlierID]
    # rmID=c(rmID,outlierID)
    # 
    # rmID=unique(rmID)
    
    ###### test directly by median
    # medL=median(adjPosL)
    # distL=adjPosL-medL
    # propL=pnorm(abs(distL), mean = fragmentInfo$fragLengthMean, sd = fragmentInfo$fragLengthSd, lower.tail = TRUE, log.p = FALSE)
    # medR=median(adjPosR)
    # distR=adjPosR-medR
    # propR=pnorm(abs(distR), mean = fragmentInfo$fragLengthMean, sd = fragmentInfo$fragLengthSd, lower.tail = TRUE, log.p = FALSE)
    # rmID=c(which(propL>0.99),which(propR>0.99)) #0.99 is too big that still allows many outliers. That might lead to decide an incorrect junction break later on
    # 
    
    ##### test the variance based on empirical variance of fragment length distribution
    medL=mean(adjPosL)
    distL=adjPosL-medL
    propL=1-sapply(distL, function(x){ min(sum(flenRef>x),sum(flenRef>x))/10000 })
    
    #propL=pnorm(abs(distL), mean = fragmentInfo$fragLengthMean, sd = fragmentInfo$fragLengthSd, lower.tail = TRUE, log.p = FALSE)
    medR=mean(adjPosR)
    distR=adjPosR-medR
    propR=1-sapply(distR, function(x){ min(sum(flenRef>x),sum(flenRef>x))/10000 })
    #propR=pnorm(abs(distR), mean = fragmentInfo$fragLengthMean, sd = fragmentInfo$fragLengthSd, lower.tail = TRUE, log.p = FALSE)
    rmID=c(which(propL<0.01),which(propR<0.01)) #0.99 is too big that still allows many outliers. That might lead to decide an incorrect junction break later on
    
    #get chrPos
    chrPosL=unlist(res$readL.chrPos)
    chrPosR=unlist(res$readR.chrPos)
    seqLenL=unlist(res$readL.seqLen)
    seqLenR=unlist(res$readR.seqLen)
    offsetL=unlist(res$readL.Pos)-unlist(res$readL.ReadPos)
    offsetR=unlist(res$readR.Pos)-unlist(res$readR.ReadPos)
    if (length(rmID)>0){
      adjPosL=adjPosL[-rmID]
      adjPosR=adjPosR[-rmID]
      chrPosL=chrPosL[-rmID]
      chrPosR=chrPosR[-rmID]
      seqLenL=seqLenL[-rmID]
      seqLenR=seqLenR[-rmID]
      offsetL=offsetL[-rmID]
      offsetR=offsetR[-rmID]
    }
    
    geneExonMatL=exonInfo[exonInfo$GENEID==gene1,]
    geneExonMatR=exonInfo[exonInfo$GENEID==gene2,]
    
    if (readStrands=="RF" || readStrands=="UN" || readStrands=="RR"){  #read1 for 3prime and read2 for 5prime
      if (length(chrPosL) > 0){
        if (junctBr$myFusionFinal$strand1[i]=="+"){
          nearest3Chrpos=min(chrPosL)
          diffL=geneExonMatL$EXONEND-nearest3Chrpos
        } else{ 
          nearest3Chrpos=max(chrPosL)
          diffL=nearest3Chrpos-geneExonMatL$EXONSTART
        }
        exIDL=which.min(abs(diffL))
        
        if (junctBr$myFusionFinal$strand2[i]=="+"){
          nearest5Chrpos=max(chrPosR)
          diffR=nearest5Chrpos-geneExonMatR$EXONSTART
        }else{
          nearest5Chrpos=min(chrPosR)
          diffR=geneExonMatR$EXONEND-nearest5Chrpos
        }
        exIDR=which.min(abs(diffR))
        diffL[exIDL]
        diffR[exIDR]
        #true breaking points
        if (junctBr$myFusionFinal$strand1[i]=="+") res$tx3br=geneExonMatL$EXONSTART[exIDL] else res$tx3br=geneExonMatL$EXONEND[exIDL]
        if (junctBr$myFusionFinal$strand2[i]=="+") res$tx5br=geneExonMatR$EXONEND[exIDR] else res$tx5br=geneExonMatR$EXONSTART[exIDR]
        
        #find the position of the true breaking points in the read distributions
        fexonL=sort(unique(c(unlist(res$readL.ExonID),geneExonMatL$EXONID[exIDL])))
        mygeneExL=geneExonMatL[match(fexonL,geneExonMatL$EXONID),]
        #chrReadPosL=unlist(res$readL.chrReadPos)
        chrReadPosL=chrPosL
        ftxReadPosL=convertChrPosGenePos(chrReadPosL,geneExonMat=mygeneExL)
        
        ftx3brPos=convertChrPosGenePos(res$tx3br,geneExonMat=mygeneExL)
        #if (junctBr$myFusionFinal$strand1[i]=="+") ftx5len=ftx5brPos-ftxReadPosL+seqLenL+offsetL else ftx5len=ftxReadPosL-ftx5brPos+seqLenL+offsetL
        if (junctBr$myFusionFinal$strand1[i]=="+") ftx3len=ftxReadPosL-ftx3brPos - offsetL + fragmentInfo$readlen else ftx3len=ftx3brPos-ftxReadPosL - offsetL + fragmentInfo$readlen
        
        fexonR=sort(unique(c(unlist(res$readR.ExonID),geneExonMatR$EXONID[exIDR])))
        mygeneExR=geneExonMatR[match(fexonR,geneExonMatR$EXONID),]
        #chrReadPosR=unlist(res$readR.chrReadPos)
        chrReadPosR=chrPosR
        ftxReadPosR=convertChrPosGenePos(chrReadPosR,geneExonMat=mygeneExR)
        ftx5brPos=convertChrPosGenePos(res$tx5br,geneExonMat=mygeneExR)
        #if (junctBr$myFusionFinal$strand2[i]=="+") ftx3len=ftxReadPosR-ftx3brPos - offsetR + fragmentInfo$readlen else ftx3len=ftx3brPos-ftxReadPosR - offsetR + fragmentInfo$readlen
        if (junctBr$myFusionFinal$strand2[i]=="+") ftx5len=ftx5brPos-ftxReadPosR+seqLenR+offsetR else ftx5len=ftxReadPosR-ftx5brPos+seqLenR+offsetR
        
        ftxLen=ftx5len+ftx3len
        
        res$fexonL=fexonL
        res$fexonR=fexonR
        
        
        res$ftxReadPosL=ftxReadPosL
        res$ftxReadPosR=ftxReadPosR
        
        res$ftx5brPos=ftx5brPos
        res$ftx5len=ftx5len
        
        res$ftx3brPos=ftx3brPos
        res$ftx3len=ftx3len
        
        res$ftxMedianLen=median(ftxLen)
        res$ftxMeanLen=mean(ftxLen)
        
        ftxLenProp=pnorm(mean(ftxLen), mean = fragmentInfo$fragLengthMean, sd = fragmentInfo$fragLengthSd, lower.tail = TRUE, log.p = FALSE)
        res$ftxLenProp=ftxLenProp
      }
      
    }else{  #read1 for 5prime and read2 for 3prime
      
      if (length(chrPosL) > 0){
        if (junctBr$myFusionFinal$strand1[i]=="+"){
          nearest5Chrpos=max(chrPosL)
          diffL=geneExonMatL$EXONEND-nearest5Chrpos
        } else{ 
          nearest5Chrpos=min(chrPosL)
          diffL=nearest5Chrpos-geneExonMatL$EXONSTART
        }
        exIDL=which.min(abs(diffL))
        
        if (junctBr$myFusionFinal$strand2[i]=="+"){
          nearest3Chrpos=min(chrPosR)
          diffR=nearest3Chrpos-geneExonMatR$EXONSTART
        }else{
          nearest3Chrpos=max(chrPosR)
          diffR=geneExonMatR$EXONEND-nearest3Chrpos
        }
        exIDR=which.min(abs(diffR))
        diffL[exIDL]
        diffR[exIDR]
        #true breaking points
        if (junctBr$myFusionFinal$strand1[i]=="+") res$tx5br=geneExonMatL$EXONEND[exIDL] else res$tx5br=geneExonMatL$EXONSTART[exIDL]
        if (junctBr$myFusionFinal$strand2[i]=="+") res$tx3br=geneExonMatR$EXONSTART[exIDR] else res$tx3br=geneExonMatR$EXONEND[exIDR]
        #find the position of the true breaking points in the read distributions
        
        fexonL=sort(unique(c(unlist(res$readL.ExonID),geneExonMatL$EXONID[exIDL])))
        mygeneExL=geneExonMatL[match(fexonL,geneExonMatL$EXONID),]
        #      chrReadPosL=unlist(res$readL.chrReadPos)
        chrReadPosL=chrPosL
        ftxReadPosL=convertChrPosGenePos(chrReadPosL,geneExonMat=mygeneExL)
        ftx5brPos=convertChrPosGenePos(res$tx5br,geneExonMat=mygeneExL)
        if (junctBr$myFusionFinal$strand1[i]=="+") ftx5len=ftx5brPos-ftxReadPosL+seqLenL+offsetL else ftx5len=ftxReadPosL-ftx5brPos+seqLenL+offsetL
        
        fexonR=sort(unique(c(unlist(res$readR.ExonID),geneExonMatR$EXONID[exIDR])))
        mygeneExR=geneExonMatR[match(fexonR,geneExonMatR$EXONID),]
        #chrReadPosR=unlist(res$readR.chrReadPos)
        chrReadPosR=chrPosR
        ftxReadPosR=convertChrPosGenePos(chrReadPosR,geneExonMat=mygeneExR)
        ftx3brPos=convertChrPosGenePos(res$tx3br,geneExonMat=mygeneExR)
        if (junctBr$myFusionFinal$strand2[i]=="+") ftx3len=ftxReadPosR-ftx3brPos - offsetR + fragmentInfo$readlen else
          ftx3len=ftx3brPos-ftxReadPosR - offsetR + fragmentInfo$readlen
        
        ftxLen=ftx5len+ftx3len
        
        res$fexonL=fexonL
        res$fexonR=fexonR
        
        res$ftxReadPosL=ftxReadPosL
        res$ftxReadPosR=ftxReadPosR
        
        res$ftx5brPos=ftx5brPos
        res$ftx5len=ftx5len
        
        res$ftx3brPos=ftx3brPos
        res$ftx3len=ftx3len
        
        res$ftxMedianLen=median(ftxLen)
        res$ftxMeanLen=mean(ftxLen)
        
        ftxLenProp=pnorm(mean(ftxLen), mean = fragmentInfo$fragLengthMean, sd = fragmentInfo$fragLengthSd, lower.tail = TRUE, log.p = FALSE)
        res$ftxLenProp=ftxLenProp
      }
    }
    res$outlierID=rmID
    res$outlierNum=length(rmID)
    junctBr$junctInfo[[i]]=res
    
    #update junctBr$myFusionFinal table
    junctBr$myFusionFinal$adjFragmentNum[i]=junctBr$myFusionFinal$fragmentNum[i]- junctBr$junctInfo[[i]]$outlierNum
    if (junctBr$myFusionFinal$adjFragmentNum[i] > 0){
      junctBr$myFusionFinal$outlierNum[i]= junctBr$junctInfo[[i]]$outlierNum
      junctBr$myFusionFinal$ftxLenProp[i]=junctBr$junctInfo[[i]]$ftxLenProp
      junctBr$myFusionFinal$ftxMedianLen[i]=junctBr$junctInfo[[i]]$ftxMedianLen
      junctBr$myFusionFinal$ftxMeanLen[i]=junctBr$junctInfo[[i]]$ftxMeanLen
      
      junctBr$myFusionFinal$ftx5LenMean[i]=mean(junctBr$junctInfo[[i]]$ftx5len)
      junctBr$myFusionFinal$ftx5LenSd[i]=sd(junctBr$junctInfo[[i]]$ftx5len)
      junctBr$myFusionFinal$ftx3LenMean[i]=mean(junctBr$junctInfo[[i]]$ftx3len)
      junctBr$myFusionFinal$ftx3LenSd[i]=sd(junctBr$junctInfo[[i]]$ftx3len)
      
      #update other information
      if (res$outlierNum > 0){
        res$crt.tx3.gap=res$tx3.gap[-rmID]
        res$crt.tx5.gap=res$tx5.gap[-rmID]
        res$crt.tx3.len=res$tx3.len[-rmID]
        res$crt.tx5.len=res$tx5.len[-rmID]
        
        if (sum(res$crt.tx3.gap==1)==0){
          res$crt.tx3.len=res$crt.tx3.len-min(res$crt.tx3.gap)+1
          res$crt.tx3.gap=res$crt.tx3.gap-min(res$crt.tx3.gap)+1
        }
        
        if (sum(res$crt.tx5.gap==1)==0){
          res$crt.tx5.len=res$crt.tx5.len-min(res$crt.tx5.gap)+1
          res$crt.tx5.gap=res$crt.tx5.gap-min(res$crt.tx5.gap)+1
        }
        
        res$crt.estFragLen=res$crt.tx3.len+res$crt.tx5.len
        
        junctBr$myFusionFinal$crt.tx5GapMean[i]=mean(res$crt.tx5.gap)
        junctBr$myFusionFinal$crt.tx3GapMean[i]=mean(res$crt.tx3.gap)
        junctBr$myFusionFinal$crt.tx5GapSd[i]=sd(res$crt.tx5.gap)
        junctBr$myFusionFinal$crt.tx3GapSd[i]=sd(res$crt.tx3.gap)
        
        junctBr$myFusionFinal$crt.tx5LenMean[i]=mean(res$crt.tx5.len)
        junctBr$myFusionFinal$crt.tx3LenMean[i]=mean(res$crt.tx3.len)
        junctBr$myFusionFinal$crt.tx5LenSd[i]=sd(res$crt.tx5.len)
        junctBr$myFusionFinal$crt.tx3LenSd[i]=sd(res$crt.tx3.len)
        
      }
    }

  }
  
  return(junctBr)
  
}




checkEndExon <-function(myFusionFinal, junctBr=NULL, anntxdb, readStrands,shrinkLen=5, type="SR"){
  ##### checking the start/ending exon of gene
  rmID=NULL
  if (type=="SR"){
    if (nrow(myFusionFinal)>0){
      exonInfo=select(anntxdb, keys=unique(c(as.character(myFusionFinal$front_gene), as.character(myFusionFinal$back_gene))), columns=c("TXNAME","EXONID","EXONSTART","EXONEND","TXSTRAND"), keytype = "GENEID")
      rmID=NULL
      for (i in 1:nrow(myFusionFinal)){
        #####5 prime
        #check exon at gene-level
        geExonMat=exonInfo[exonInfo$GENEID==myFusionFinal$gene1[i],]
        exStartID=which(geExonMat$EXONSTART==min(geExonMat$EXONSTART))
        exEndID=which(geExonMat$EXONEND==max(geExonMat$EXONEND))
        diffS=unlist(lapply(myFusionFinal$brchposEx5[i],function(x) x-geExonMat$EXONSTART))
        diffE=unlist(lapply(myFusionFinal$brchposEx5[i],function(x) x-geExonMat$EXONEND))
        if (geExonMat$TXSTRAND[1]=="-") diff=diffS  else diff=diffE
        if (min(abs(diff)) < shrinkLen) {
          myExID=which(abs(diff) < shrinkLen)
          if (geExonMat$TXSTRAND[1]=="-") exBoundID=exStartID else exBoundID=exEndID
          if (sum(is.na(match(myExID,exBoundID)))==0){
            rmID=c(rmID,i)
          }
        } #else{} # if so, the exon of the junction can be detected as the follow: choose the other diff(diffS/diffE) and select the start/end of the next/previous exon of the transcripts in the genes. However, that mean this one is not the start/end exon of the gene. Thus we don't need to determine the exon ID anymore
        ####3 prime
        #check exon at gene-level
        geExonMat=exonInfo[exonInfo$GENEID==myFusionFinal$gene2[i],]
        exStartID=which(geExonMat$EXONSTART==min(geExonMat$EXONSTART))
        exEndID=which(geExonMat$EXONEND==max(geExonMat$EXONEND))
        diffS=unlist(lapply(myFusionFinal$brchposEx3[i],function(x) x-geExonMat$EXONSTART))
        diffE=unlist(lapply(myFusionFinal$brchposEx3[i],function(x) x-geExonMat$EXONEND))
        if (geExonMat$TXSTRAND[1]=="+") diff=diffS else diff=diffE
        
        if (min(abs(diff)) < shrinkLen) {
          myExID=which(abs(diff) < shrinkLen)
          if (geExonMat$TXSTRAND[1]=="+") exBoundID=exStartID else exBoundID=exEndID
          if (sum(is.na(match(myExID,exBoundID)))==0){
            rmID=c(rmID,i)
          }
        }#else{} # if so, the exon of the junction can be detected as the follow: choose the other diff(diffE/diffS) and select the start/end of the next/previous exon of the transcripts in the genes. However, that mean this one is not the start/end exon of the gene. Thus we don't need to determine the exon ID anymore
      }
    }
  }
  
  if (type=="MR"){
    if (nrow(myFusionFinal)>0){
      exonInfo=select(anntxdb, keys=unique(c(as.character(myFusionFinal$gene1), as.character(myFusionFinal$gene2))), columns=c("TXNAME","EXONID","EXONSTART","EXONEND","TXSTRAND"), keytype = "GENEID")
      
      for (i in 1:nrow(myFusionFinal)){
        myfge=myFusionFinal$fusionName[i]
        myjbr=junctBr$junctInfo[[myfge]]
        ##### check 5 prime
        if (readStrands=="RF" || readStrands=="UN" || readStrands=="RR"){  #read1 for 3prime and read2 for 5prime
          geExonMat=exonInfo[exonInfo$GENEID==myFusionFinal$gene2[i],]
        }else{ #read1 for 5prime and read2 for 3prime
          geExonMat=exonInfo[exonInfo$GENEID==myFusionFinal$gene1[i],]
        }
        
        exStartID=which(geExonMat$EXONSTART==min(geExonMat$EXONSTART))
        exEndID=which(geExonMat$EXONEND==max(geExonMat$EXONEND))
        
        if (geExonMat$TXSTRAND[1]=="-") exBoundID=exStartID else exBoundID=exEndID
        if (geExonMat$TXSTRAND[1]=="-") myExID=geExonMat$EXONID[which(geExonMat$EXONSTART==myjbr$tx5br)] else myExID=geExonMat$EXONID[which(geExonMat$EXONEND==myjbr$tx5br)]
        
        if(sum(!is.na(match(unique(geExonMat$EXONID[exBoundID]),unique(myExID))))>0) rmID=c(rmID,i)
        
        
        ##### check 3 prime
        if (readStrands=="RF" || readStrands=="UN" || readStrands=="RR"){  #read1 for 3prime and read2 for 5prime
          geExonMat=exonInfo[exonInfo$GENEID==myFusionFinal$gene1[i],]
          
        }else{ #read1 for 5prime and read2 for 3prime
          geExonMat=exonInfo[exonInfo$GENEID==myFusionFinal$gene2[i],]
        }
        
        exStartID=which(geExonMat$EXONSTART==min(geExonMat$EXONSTART))
        exEndID=which(geExonMat$EXONEND==max(geExonMat$EXONEND))
        if (geExonMat$TXSTRAND[1]=="+") exBoundID=exStartID else exBoundID=exEndID
        if (geExonMat$TXSTRAND[1]=="+") myExID=geExonMat$EXONID[which(geExonMat$EXONSTART==myjbr$tx3br)] else myExID=geExonMat$EXONID[which(geExonMat$EXONEND==myjbr$tx3br)]
        if(sum(!is.na(match(unique(geExonMat$EXONID[exBoundID]),unique(myExID))))>0) rmID=c(rmID,i)
        
      }
    }
  }
  
  return(rmID)
}








