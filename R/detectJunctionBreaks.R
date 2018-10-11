############################################################
##### detect junction break positions from mapped reads
detectJunctionBreaks <-function(fgeList,inPath,feq, feqFgeMap, anntxdb, readStrands="UN", shrinkLen=3){
  myFusionFinal=fgeList
  ##### input fusion reads
  #load fusion reads
  cat("\n Get mapped fusion reads")
  frfiles=list.files(inPath,paste(readStrands,"_fusionMappedReadsChunk_*",sep=""))
  fusionRead=NULL;
  fsizes=NULL;
  for (i in 1:length(frfiles)){
    tmpDat=read.csv(paste(inPath,"/",frfiles[i],sep=""), header =FALSE, sep="\t")
    fusionRead=rbind(fusionRead,tmpDat)
    fsizes=c(fsizes,nrow(tmpDat))
  }
  fsizeLadder=cumsum(fsizes)
  colnames(fusionRead)=c("read1","read2","read1Pos","read2Pos","seq1Pos","seq2Pos","seq1Len","seq2Len")
  #dim(fusionRead)
  cat("\n The number of mapped fusion reads: ",nrow(fusionRead))
  fre1=as.character(fusionRead$read1)
  fre1=trimws(fre1)
  fre2=as.character(fusionRead$read2)
  fre2=trimws(fre2)
  read1Pos=as.character(fusionRead$read1Pos)
  read2Pos=as.character(fusionRead$read2Pos)
  seq1Pos=as.character(fusionRead$seq1Pos)
  seq2Pos=as.character(fusionRead$seq2Pos)
  seq1Len=as.character(fusionRead$seq1Len)
  seq2Len=as.character(fusionRead$seq2Len)
  

    #load fragment info - this is not neccessary thi moment
  fragmentInfo=read.csv(paste(inPath,"/fragmentInfo.txt",sep=""), header =TRUE, sep="\t")
  readLen=fragmentInfo[1,1]
  #load the feq file
  feqRaw=read.csv(paste(inPath,"/feq_",readStrands,".txt",sep=""), header =TRUE, sep="\t")
  
  # #Keep only reads and feqs relating to myFusionFinal
  # allTx=select(anntxdb, keys=unique(c(as.character(myFusionFinal$gene5p),as.character(myFusionFinal$gene3p))), columns=c("TXNAME"), keytype = "GENEID")
  # feqRaw=feqRaw[as.character(feqRaw$Transcript) %in% as.character(allTx$TXNAME),]
  
  feqRead1=feqRaw[feqRaw$Read==1,]
  feqRead2=feqRaw[feqRaw$Read==2,]
  #get feq-fge map
  # feqFgeMap=feqInfo$feqFgeMap
  # feq=feqInfo$feq
  
  
  
  
  cat("\n Preparing other information ...")
  feqFtxMap1=tapply(as.character(feqRead1$Transcript),feqRead1$Feq,c)
  feqRead1Name=unlist(lapply(feqFtxMap1,function(x) paste(x,collapse =" ")))
  feqFtxMap2=tapply(as.character(feqRead2$Transcript),feqRead2$Feq,c)
  feqRead2Name=unlist(lapply(feqFtxMap2,function(x) paste(x,collapse =" ")))
  
  read1Pos=lapply(read1Pos,function(x) as.integer(unlist(strsplit(x," "))))
  read2Pos=lapply(read2Pos,function(x) as.integer(unlist(strsplit(x," "))))
  seq1Pos=lapply(seq1Pos,function(x) as.integer(unlist(strsplit(x," "))))
  seq2Pos=lapply(seq2Pos,function(x) as.integer(unlist(strsplit(x," "))))
  seq1Len=lapply(seq1Len,function(x) as.integer(unlist(strsplit(x," "))))
  seq2Len=lapply(seq2Len,function(x) as.integer(unlist(strsplit(x," "))))
  
  
  #get some annotations
  exonInfo=select(anntxdb, keys=unique(c(as.character(myFusionFinal$gene1),as.character(myFusionFinal$gene2))), columns=c("TXNAME","EXONID","EXONSTART","EXONEND","TXSTRAND"), keytype = "GENEID")
  txToGene=select(anntxdb, keys=unique(as.character(feqRaw[,1])), columns=c("GENEID","TXCHROM"), keytype = "TXNAME")


  cat("\n Detect positions of junction breaks")
  brpos5.start=brpos5.end=brpos3.start=brpos3.end=NULL
  genebrpos5.start=genebrpos5.end=genebrpos3.start=genebrpos3.end=NULL
  fragmentMean=fragmentSd=fragmentNum=fragmentTest=NULL;
  tx5LenMean=tx3LenMean=tx5LenSd=tx3LenSd=NULL;
  tx5GapMean=tx3GapMean=tx5GapSd=tx3GapSd=NULL;
  flen5=flen3=NULL;
  splitRNum5=splitRNum3=exNum1=exNum2=NULL;
  nondupCount=NULL;
  
  junctInfo=list();
  
  for (i in 1:nrow(myFusionFinal)){
    fgeName=as.character(myFusionFinal$name12)[i]

    gene1=as.character(myFusionFinal$gene1[i])
    gene2=as.character(myFusionFinal$gene2[i])
    myfeqID=feqFgeMap[[fgeName]]

    readL.seqLen=readR.seqLen=readID=readL.chrPos=readR.chrPos=readL.Pos=readR.Pos=readL.GenePos=readR.GenePos=readL.tx=readR.tx=list();
    readL.chrReadPos=readR.chrReadPos=readR.ReadPos=readL.ReadPos=list();
    readL.ExonID=readR.ExonID=list();

    for (j in 1:length(myfeqID)){
      feqName=names(feq[myfeqID[j]])
      #keepID=which(fre1==feqRead1Name[myfeqID[j]] & fre2==feqRead2Name[myfeqID[j]])
      keepID1=which(fre1==feqRead1Name[myfeqID[j]])
      keepID=keepID1[which(fre2[keepID1]==feqRead2Name[myfeqID[j]])]

      readID[[feqName]]=keepID
      
      xl=unlist(feqFtxMap1[myfeqID[j]])
      xr=unlist(feqFtxMap2[myfeqID[j]])
      

      IDl=which(txToGene[match(xl,txToGene$TXNAME),]$GENEID==gene1)
      IDr=which(txToGene[match(xr,txToGene$TXNAME),]$GENEID==gene2)
      #select only one because the others will locate to the same posistion
      posl=IDl[1]
      posr=IDr[1]

      readposl=sapply(read1Pos[keepID],function(x) x[posl])
      readposr=sapply(read2Pos[keepID],function(x) x[posr])
      txposl=sapply(seq1Pos[keepID],function(x) x[posl])
      txposr=sapply(seq2Pos[keepID],function(x) x[posr])

      txSeqlenl=sapply(seq1Len[keepID],function(x) x[posl])
      txSeqlenr=sapply(seq2Len[keepID],function(x) x[posr])

      #reduce length of match kmer shrinkLen=3 or prop=0.25^3=0.015625
      readposl=readposl+shrinkLen
      readposr=readposr+shrinkLen
      txposl=txposl+shrinkLen
      txposr=txposr+shrinkLen
      txSeqlenl=txSeqlenl-shrinkLen
      txSeqlenr=txSeqlenr-shrinkLen
      
      #save the information
      readL.Pos[[feqName]]=txposl
      readR.Pos[[feqName]]=txposr
      readL.ReadPos[[feqName]]=readposl
      readR.ReadPos[[feqName]]=readposr
      readL.seqLen[[feqName]]=txSeqlenl
      readR.seqLen[[feqName]]=txSeqlenr
      

      
      #find position at chromosome level  
      txnamel=xl[posl]
      txchrposl=convertChrPos(txnamel,txposl,exonInfo)
      readL.chrPos[[feqName]]=txchrposl
      txnamer=xr[posr]
      txchrposr=convertChrPos(txnamer,txposr,exonInfo)
      readR.chrPos[[feqName]]=txchrposr
      
      readL.chrReadPos[[feqName]]=convertChrPos(txnamel,readposl,exonInfo)
      readR.chrReadPos[[feqName]]=convertChrPos(txnamer,readposr,exonInfo)
      
      
      #find the position at gene-level
      readL.GenePos[[feqName]]=convertChrPosGenePos(txchrposl,gene1,exonInfo)
      readR.GenePos[[feqName]]=convertChrPosGenePos(txchrposr,gene2,exonInfo)

      readL.ExonID[[feqName]]=detectExonID3(txnamel,txposl,txSeqlenl,exonInfo)
      readR.ExonID[[feqName]]=detectExonID3(txnamer,txposr,txSeqlenr, exonInfo)

      readL.tx[[feqName]]=txnamel
      readR.tx[[feqName]]=txnamer
      

    }
    
    res=list(readL.ReadPos=readL.ReadPos,readR.ReadPos=readR.ReadPos,readL.seqLen=readL.seqLen, readR.seqLen=readR.seqLen,readL.chrPos=readL.chrPos,readR.chrPos=readR.chrPos,readL.Pos=readL.Pos,readR.Pos=readR.Pos, readL.GenePos=readL.GenePos,readR.GenePos=readR.GenePos,readL.ExonID=readL.ExonID,readR.ExonID=readR.ExonID,readID=readID,readL.tx=readL.tx,readR.tx=readR.tx,readL.chrReadPos=readL.chrReadPos,readR.chrReadPos=readR.chrReadPos)

    res$new.readL.GenePos=res$new.readR.GenePos=list()
    mygeneExL=exonInfo[match(sort(unique(unlist(res$readL.ExonID))),exonInfo$EXONID),]
    for (j in 1:length(myfeqID)) res$new.readL.GenePos[[names(res$readL.chrPos[j])]]=convertChrPosGenePos(res$readL.chrPos[[j]],geneExonMat=mygeneExL)
    mygeneExR=exonInfo[match(sort(unique(unlist(res$readR.ExonID))),exonInfo$EXONID),]
    for (j in 1:length(myfeqID)) res$new.readR.GenePos[[names(res$readR.chrPos[j])]]=convertChrPosGenePos(res$readR.chrPos[[j]],geneExonMat=mygeneExR)
    
    ######   
    exNum1=c(exNum1,length(unique(unlist(res$readL.ExonID))))
    exNum2=c(exNum2,length(unique(unlist(res$readR.ExonID))))
    #update trueStartPos of fusion-gene
    nondupCount=c(nondupCount,length(unique(paste(unlist(res$readL.chrReadPos),unlist(res$readR.chrReadPos),sep=""))))
    readDirectTag=c(0,0)
    
    res$adj.readL.GenePos=res$new.readL.GenePos
    res$adj.readR.GenePos=res$new.readR.GenePos
    
    ### In real data, we always work on RF direction even for UN (unstranded), so 5prime is in the Right side (read2) and 3prime is in the leftSide (read1)
    # If RR happens, it might be a reorder of genes or due to protocols that should be investigate further.
    if (readStrands=="RF" || readStrands=="UN" || readStrands=="RR"){ #read1 for 3prime and read2 for 5prime
      readDirectTag=c(3,5)
      flen3=c(flen3,getGeneLen(geneExonMat=mygeneExL))
      flen5=c(flen5,getGeneLen(geneExonMat=mygeneExR))
      if (myFusionFinal$strand1[i] == "+"){
        brpos3.start=c(brpos3.start,min(unlist(res$readL.chrPos)))
        brpos3.end=c(brpos3.end,max(unlist(res$readL.chrPos)))
        genebrpos3.start=c(genebrpos3.start,min(unlist(res$adj.readL.GenePos)))
        genebrpos3.end=c(genebrpos3.end,max(unlist(res$adj.readL.GenePos)))
        mySR3=0
        for (myid in 1:length(res$adj.readL.GenePos)) mySR3=mySR3+sum(res$adj.readL.GenePos[[myid]]-(res$readL.Pos[[myid]]-res$readL.ReadPos[[myid]])<genebrpos3.start[length(genebrpos3.start)]-shrinkLen)
        splitRNum3=c(splitRNum3,mySR3)
      } else {
        brpos3.start=c(brpos3.start,max(unlist(res$readL.chrPos)))
        brpos3.end=c(brpos3.end,min(unlist(res$readL.chrPos)))
        genebrpos3.start=c(genebrpos3.start,max(unlist(res$adj.readL.GenePos)))
        genebrpos3.end=c(genebrpos3.end,min(unlist(res$adj.readL.GenePos)))
        mySR3=0
        for (myid in 1:length(res$adj.readL.GenePos)) mySR3=mySR3+sum(res$adj.readL.GenePos[[myid]]+(res$readL.Pos[[myid]]-res$readL.ReadPos[[myid]])>genebrpos3.start[length(genebrpos3.start)]+shrinkLen)
        splitRNum3=c(splitRNum3,mySR3)
      }
      
      if (myFusionFinal$strand2[i] == "+"){
        for (myid in 1:length(res$adj.readR.GenePos)) res$adj.readR.GenePos[[myid]]=res$adj.readR.GenePos[[myid]]+res$readR.seqLen[[myid]]-1
        brpos5.start=c(brpos5.start,max(unlist(res$readR.chrPos)))
        brpos5.end=c(brpos5.end,min(unlist(res$readR.chrPos)))
        genebrpos5.start=c(genebrpos5.start,max(unlist(res$adj.readR.GenePos)))
        genebrpos5.end=c(genebrpos5.end,min(unlist(res$adj.readR.GenePos)))
        mySR5=0
        for (myid in 1:length(res$adj.readR.GenePos)) mySR5=mySR5+sum(res$adj.readR.GenePos[[myid]]-(res$readR.Pos[[myid]]-res$readR.ReadPos[[myid]])-(res$readR.seqLen[[myid]]-1)+fragmentInfo$readlen-1>genebrpos5.start[length(genebrpos5.start)]+shrinkLen)
        splitRNum5=c(splitRNum5,mySR5)
      } else{
        for (myid in 1:length(res$adj.readR.GenePos)) res$adj.readR.GenePos[[myid]]=res$adj.readR.GenePos[[myid]]-res$readR.seqLen[[myid]]+1
        brpos5.start=c(brpos5.start,min(unlist(res$readR.chrPos)))
        brpos5.end=c(brpos5.end,max(unlist(res$readR.chrPos)))
        genebrpos5.start=c(genebrpos5.start,min(unlist(res$adj.readR.GenePos)))
        genebrpos5.end=c(genebrpos5.end,max(unlist(res$adj.readR.GenePos)))
        mySR5=0
        for (myid in 1:length(res$adj.readR.GenePos)) mySR5=mySR5+sum(res$adj.readR.GenePos[[myid]]+(res$readR.Pos[[myid]]-res$readR.ReadPos[[myid]])+(res$readR.seqLen[[myid]]-1)-fragmentInfo$readlen+1<genebrpos5.start[length(genebrpos5.start)]-shrinkLen)
        splitRNum5=c(splitRNum5,mySR5)
        
      }
    } else{
      readDirectTag=c(5,3)
      flen5=c(flen5,getGeneLen(geneExonMat=mygeneExL))
      flen3=c(flen3,getGeneLen(geneExonMat=mygeneExR))
      #Process for FR, FF: read1 for 5prime and read2 for 3prime. If FF happens, it might be a reorder of genes or due to protocols that should be investigate further.
      if (myFusionFinal$strand1[i] == "+"){
        for (myid in 1:length(res$adj.readL.GenePos)) res$adj.readL.GenePos[[myid]]=res$adj.readL.GenePos[[myid]]+res$readL.seqLen[[myid]]-1
        brpos5.start=c(brpos5.start,max(unlist(res$readL.chrPos)))
        brpos5.end=c(brpos5.end,min(unlist(res$readL.chrPos)))
        genebrpos5.start=c(genebrpos5.start,max(unlist(res$adj.readL.GenePos)))
        genebrpos5.end=c(genebrpos5.end,min(unlist(res$adj.readL.GenePos)))
        mySR5=0
        for (myid in 1:length(res$adj.readL.GenePos)) mySR5=mySR5+sum(res$adj.readL.GenePos[[myid]]-(res$readL.Pos[[myid]]-res$readL.ReadPos[[myid]])-(res$readL.seqLen[[myid]]-1)+fragmentInfo$readlen-1>genebrpos5.start[length(genebrpos5.start)]+shrinkLen)
        splitRNum5=c(splitRNum5,mySR5)
      }else{
        for (myid in 1:length(res$adj.readL.GenePos)) res$adj.readL.GenePos[[myid]]=res$adj.readL.GenePos[[myid]]-res$readL.seqLen[[myid]]+1
        brpos5.start=c(brpos5.start,min(unlist(res$readL.chrPos)))
        brpos5.end=c(brpos5.end,max(unlist(res$readL.chrPos)))
        genebrpos5.start=c(genebrpos5.start,min(unlist(res$adj.readL.GenePos)))
        genebrpos5.end=c(genebrpos5.end,max(unlist(res$adj.readL.GenePos)))
        mySR5=0
        for (myid in 1:length(res$adj.readL.GenePos)) mySR5=mySR5+sum(res$adj.readL.GenePos[[myid]]+(res$readL.Pos[[myid]]-res$readL.ReadPos[[myid]])+(res$readL.seqLen[[myid]]-1)-fragmentInfo$readlen+1<genebrpos5.start[length(genebrpos5.start)]-shrinkLen)
        splitRNum5=c(splitRNum5,mySR5)
      }
      
      if (myFusionFinal$strand2[i] == "+"){
        brpos3.start=c(brpos3.start,min(unlist(res$readR.chrPos)))
        brpos3.end=c(brpos3.end,max(unlist(res$readR.chrPos)))
        genebrpos3.start=c(genebrpos3.start,min(unlist(res$adj.readR.GenePos)))
        genebrpos3.end=c(genebrpos3.end,max(unlist(res$adj.readR.GenePos)))
        mySR3=0
        for (myid in 1:length(res$adj.readR.GenePos)) mySR3=mySR3+sum(res$adj.readR.GenePos[[myid]]-(res$readR.Pos[[myid]]-res$readR.ReadPos[[myid]])<genebrpos3.start[length(genebrpos3.start)]-shrinkLen)
        splitRNum3=c(splitRNum3,mySR3)
      } else{
        brpos3.start=c(brpos3.start,max(unlist(res$readR.chrPos)))
        brpos3.end=c(brpos3.end,min(unlist(res$readR.chrPos)))
        genebrpos3.start=c(genebrpos3.start,max(unlist(res$adj.readR.GenePos)))
        genebrpos3.end=c(genebrpos3.end,min(unlist(res$adj.readR.GenePos)))
        mySR3=0
        for (myid in 1:length(res$adj.readR.GenePos)) mySR3=mySR3+sum(res$adj.readR.GenePos[[myid]]+(res$readR.Pos[[myid]]-res$readR.ReadPos[[myid]])>genebrpos3.start[length(genebrpos3.start)]+shrinkLen)
        splitRNum3=c(splitRNum3,mySR3)
      }
    }
    
      #extract fragment distributions
    if (readDirectTag[1]==3){
      tx3.gap=abs(unlist(res$adj.readL.GenePos)-genebrpos3.start[length(genebrpos3.start)])+1
      tx5.gap=abs(unlist(res$adj.readR.GenePos)-genebrpos5.start[length(genebrpos5.start)])+1
      tx3.len=tx3.gap+unlist(res$readL.seqLen)
      tx5.len=tx5.gap+unlist(res$readR.seqLen)
      #fragLen=tx3.len+tx5.len+unlist(res$readL.seqLen)+unlist(res$readR.seqLen)
    }else{
      tx5.gap=abs(unlist(res$adj.readL.GenePos)-genebrpos5.start[length(genebrpos5.start)])+1
      tx3.gap=abs(unlist(res$adj.readR.GenePos)-genebrpos3.start[length(genebrpos3.start)])+1
      
      tx5.len=tx5.gap+unlist(res$readL.seqLen)
      tx3.len=tx3.gap+unlist(res$readR.seqLen)
      #fragLen=tx3.len+tx5.len+unlist(res$readR.seqLen)+unlist(res$readL.seqLen)
    }
    #get fragment length
    #plot(tx3.len,tx5.len)

    res$tx3.gap=tx3.gap
    res$tx5.gap=tx5.gap
    res$tx3.len=tx3.len
    res$tx5.len=tx5.len
    
    res$estFragLen=tx3.len+tx5.len
    
    tx5GapMean=c(tx5GapMean,mean(tx5.gap))
    tx5GapSd=c(tx5GapSd,sd(tx5.gap))
    tx3GapMean=c(tx3GapMean,mean(tx3.gap))
    tx3GapSd=c(tx3GapSd,sd(tx3.gap))

    tx5LenMean=c(tx5LenMean,mean(tx5.len))
    tx5LenSd=c(tx5LenSd,sd(tx5.len))
    tx3LenMean=c(tx3LenMean,mean(tx3.len))
    tx3LenSd=c(tx3LenSd,sd(tx3.len))

    fragmentMean=c(fragmentMean,mean(res$estFragLen))
    fragmentSd=c(fragmentSd,sd(res$estFragLen))
    fragmentNum=c(fragmentNum,length(res$estFragLen))
    # do a statistic test
    myttest=pnorm(fragmentMean[length(fragmentMean)],mean=fragmentInfo$fragLengthMean,sd=fragmentInfo$fragLengthSd, lower.tail = TRUE)
    fragmentTest=c(fragmentTest,myttest)
 
    junctInfo[[myFusionFinal$fusionName[i]]]=res;
  }  
  

  ###### now we work on gene5 and gene3 only
  myFusionFinal$nondupCount=nondupCount
  myFusionFinal$exNum1=exNum1
  myFusionFinal$exNum2=exNum2
  myFusionFinal$fragmentMean=fragmentMean
  myFusionFinal$fragmentSd=fragmentSd
  myFusionFinal$fragmentNum=fragmentNum
  myFusionFinal$fragmentTest=fragmentTest
  myFusionFinal$tx5GapMean=tx5GapMean
  myFusionFinal$tx5GapSd=tx5GapSd
  myFusionFinal$tx3GapMean=tx3GapMean
  myFusionFinal$tx3GapSd=tx3GapSd

  myFusionFinal$tx5LenMean=tx5LenMean
  myFusionFinal$tx5LenSd=tx5LenSd
  myFusionFinal$tx3LenMean=tx3LenMean
  myFusionFinal$tx3LenSd=tx3LenSd
  
    
  myFusionFinal$flen5=flen5
  myFusionFinal$flen3=flen3
  myFusionFinal$splitRNum5=splitRNum5
  myFusionFinal$splitRNum3=splitRNum3
  
  
  myFusionFinal$brpos5.start=brpos5.start
  myFusionFinal$brpos3.start=brpos3.start
  myFusionFinal$brpos5.end=brpos5.end
  myFusionFinal$brpos3.end=brpos3.end
  
  myFusionFinal$genebrpos5.start=genebrpos5.start
  myFusionFinal$genebrpos5.end=genebrpos5.end
  myFusionFinal$genebrpos3.start=genebrpos3.start
  myFusionFinal$genebrpos3.end=genebrpos3.end
  
  myFusionFinal$genebrpos5.rg=abs(myFusionFinal$genebrpos5.start-myFusionFinal$genebrpos5.end)+1
  myFusionFinal$genebrpos3.rg=abs(myFusionFinal$genebrpos3.start-myFusionFinal$genebrpos3.end)+1
  
  junctDist.mat=cbind(abs(myFusionFinal$brpos5.start-myFusionFinal$brpos3.start),
            abs(myFusionFinal$brpos5.start-myFusionFinal$brpos3.end),
            abs(myFusionFinal$brpos5.end-myFusionFinal$brpos3.start),
            abs(myFusionFinal$brpos5.end-myFusionFinal$brpos3.end))
  junctDist=unlist(apply(junctDist.mat,1,min))
  myFusionFinal$junctDist=junctDist
  
 
  cat("\n Get reads supporting constituent genes")
  txeq=read.csv(paste(inPath,"/rawCount.txt",sep=""), header =TRUE, sep="\t")
  #txeq$geneID=geneAnno[match(as.character(txeq$Transcript),geneAnno[,1]),6] # risky: do not pass geneAnno but still use here
  txeq$geneID=txToGene$GENEID[match(as.character(txeq$Transcript),as.character(txToGene$TXNAME))]
    
  geeq=txeq[,c(3,7,8)]
  geeq=geeq[!duplicated(geeq),]
  geCount=tapply(geeq$Count,geeq$geneID,sum)
  myFusionFinal$gene1Count=geCount[match(as.character(myFusionFinal$gene1),names(geCount))]
  myFusionFinal$gene2Count=geCount[match(as.character(myFusionFinal$gene2),names(geCount))]
  myFusionFinal$gene1Count[is.na(myFusionFinal$gene1Count)]=0
  myFusionFinal$gene2Count[is.na(myFusionFinal$gene2Count)]=0
  
  #number of feq supporting a fge
  #feqNum=sapply(as.character(myFusionFinal$name12),function(mykey) length(feqFgeMap[[mykey]]))
  matchID=match(as.character(myFusionFinal$name12),names(feqFgeMap))
  myfeqID=feqFgeMap[matchID]
  feqNum=sapply(myfeqID, length)
  
  myFusionFinal$feqNum=feqNum

  return(list(myFusionFinal=myFusionFinal,junctInfo=junctInfo,fsizeLadder=fsizeLadder))
  }
  
  