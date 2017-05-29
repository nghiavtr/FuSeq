############################################################
#####post process split reads
############################################################

postProcessSplitRead <-function(inPath, anntxdb, FuSeq.SR, FuSeq.MR, txFastaFile, FuSeq.params, shrinkLen=5){
cat("\n Post processing split reads (SR)...")
##### post processing with information from split reads
myFusionTmp=FuSeq.SR$myFusionFinal
fragmentInfo=FuSeq.SR$fragmentInfo
fragDist=FuSeq.SR$fragDist


mappedFge=FuSeq.MR$feqInfo$fgeList
readStrands=FuSeq.params$readStrands
if (readStrands=="RF" || readStrands=="UN" || readStrands=="RR"){
  mappedFge$gene5p=mappedFge$gene2
  mappedFge$gene3p=mappedFge$gene1
  mappedFge$chrom5p=mappedFge$chrom2
  mappedFge$chrom3p=mappedFge$chrom1
  mappedFge$strand5p=mappedFge$strand2
  mappedFge$strand3p=mappedFge$strand1

}else{
  mappedFge$gene5p=mappedFge$gene1
  mappedFge$gene3p=mappedFge$gene2
  mappedFge$chrom5p=mappedFge$chrom1
  mappedFge$chrom3p=mappedFge$chrom2
  mappedFge$strand5p=mappedFge$strand1
  mappedFge$strand3p=mappedFge$strand2
}
mappedFge$fusionName=paste(mappedFge$gene5p,mappedFge$gene3p,sep="-")


myFusionTmp2=myFusionTmp
#If we limit only splicing sites at exon boundary
myFusionTmp2=myFusionTmp2[myFusionTmp2$ssEnd<=2,]
myFusionTmp2=myFusionTmp2[myFusionTmp2$ssStart<=2,]
matchID=match(as.character(myFusionTmp2$name12),as.character(mappedFge$fusionName))

myFusionTmp2$total21=myFusionTmp2$tx21Count
myFusionTmp2$total21[which(!is.na(matchID))]=myFusionTmp2$total21[which(!is.na(matchID))]+mappedFge$name21Count[na.omit(matchID)]
myFusionTmp2$totalCount=myFusionTmp2$adjtx12Count
myFusionTmp2$totalCount[which(!is.na(matchID))]=myFusionTmp2$totalCount[which(!is.na(matchID))]+mappedFge$correctedCount[na.omit(matchID)]
myFusionTmp2$mappedCrtCount=myFusionTmp2$totalCount-myFusionTmp2$adjtx12Count
myFusionTmp2$mappedCount=myFusionTmp2$mappedCrtCount
myFusionTmp2$mappedCount=0
myFusionTmp2$mappedCount[which(!is.na(matchID))]=mappedFge$supportCount[na.omit(matchID)]

#filter again by inverted direction fusion genes
myFusionTmp2=myFusionTmp2[myFusionTmp2$total21 <= myFusionTmp2$totalCount*0.01,]


##### process cases: geneA-geneB vs geneA-geneC where geneB and geneC are paralogs or overlapping
myFusion=myFusionTmp2
#do some filters here
# totalCount here is the sum of adjusted counts from mapped reads and split reads
myFusion=myFusion[myFusion$totalCount>=1,]

myFusion$note3=myFusion$note5=""
myDup=duplicated(myFusion$name12)
myFusionNoDup=myFusion[!myDup,]
dupge1=table(myFusionNoDup$gene1)
dupge2=table(myFusionNoDup$gene2)
myFusion=cbind(myFusion,as.integer(dupge1[match(as.character(myFusion$gene1),names(dupge1))]))
myFusion=cbind(myFusion,as.integer(dupge2[match(as.character(myFusion$gene2),names(dupge2))]))
colnames(myFusion)[c(ncol(myFusion)-1,ncol(myFusion))]=c("dupge1f","dupge2f")

exonInfo=select(anntxdb, keys=unique(c(as.character(myFusion$gene1),as.character(myFusion$gene2))), columns=c("TXNAME","EXONID","EXONSTART","EXONEND","TXSTRAND"), keytype = "GENEID")


frontDupge=unique(as.character(myFusion$gene1[which(myFusion$dupge1f>1)]))
rmFusion=NULL;
if (length(frontDupge)>0)
  for (myge in frontDupge){
    dupFusion=myFusion[which(as.character(myFusion$gene1)==myge),]
    mygene2=unique(as.character(dupFusion$gene2))
    #check if two adjacent genes are overlapped
    myEx=exonInfo[which(!is.na(match(exonInfo$GENEID,mygene2))),]
    myEx=myEx[order(myEx$EXONSTART, decreasing = FALSE),] #keep genes in the order of exonstart
    mygene2=unique(as.character(myEx$GENEID))
    brStart=tapply(myEx$EXONSTART,myEx$GENEID,min)
    brEnd=tapply(myEx$EXONEND,myEx$GENEID,max)
    
    olid1=olid2=NULL
    for (i in 1:(length(brStart)-1)){
      for (j in (i+1):length(brStart)){
        isParalog=FALSE
        myStatusij=(brStart[i]-brStart[j])/1e8*(brEnd[i]-brStart[j])
        myStatusji=(brStart[j]-brStart[i])/1e8*(brEnd[j]-brStart[i])
        if (myStatusij < 0 | myStatusji <0) isParalog=TRUE
        if (!isParalog){
          par1=c(as.character(mygene2[i]),geneParalog[which(geneParalog[,1]==as.character(mygene2[i])),2])
          par2=c(as.character(mygene2[j]),geneParalog[which(geneParalog[,1]==as.character(mygene2[j])),2])
          if (length(intersect(par1,par2))>0) isParalog=TRUE
        }
        if (isParalog){
          olid1=c(olid1,i)
          olid2=c(olid2,j)
        }
      }
    }
    if (length(olid1)>0){
      olGroup=seq(brStart)
      olGroup[olid2]=olid1
      olGroupU=rep(-1,length(olGroup))
      while(sum(olGroup!=olGroupU)>0){
        olGroupU=olGroup
        olGroup=olGroup[olGroup]
      }
      olGroupID=unique(olGroup)
      for (myol in olGroupID){#select only one
        keepID=which(olGroup==myol)
        if (length(keepID) >1){
          olFusionName=paste(myge,mygene2[keepID],sep="-")
          olFusion=dupFusion[match(olFusionName,dupFusion$name12),]
          selectID=which.max(olFusion$totalCount)
          #check if no fusion gene has higher total count
          if (sum(olFusion$totalCount[selectID] > olFusion$totalCount[-selectID])==0){
            #so we need another way to select the best candidate, we select the one closer geneDist
            selectID=which.min(olFusion$geneDist)
          }
          #update to the removal list
          keepID=which(!is.na(match(myFusion$name12,olFusionName)))
          myFusion$note3[keepID]=paste(myFusion$note3[keepID],paste(olFusionName,collapse = ","),sep=";")
        }        
      }
    }
  }

backDupge=unique(as.character(myFusion$gene2[which(myFusion$dupge2f>1)]))
rmFusion=NULL;
if (length(backDupge)>0)
  for (myge in backDupge){
    dupFusion=myFusion[which(as.character(myFusion$gene2)==myge),]
    mygene1=unique(as.character(dupFusion$gene1))
    #check if two adjacent genes are overlapped
    myEx=exonInfo[which(!is.na(match(exonInfo$GENEID,mygene1))),]
    myEx=myEx[order(myEx$EXONSTART, decreasing = FALSE),] #keep genes in the order of exonstart
    mygene1=unique(as.character(myEx$GENEID))
    brStart=tapply(myEx$EXONSTART,myEx$GENEID,min)
    brEnd=tapply(myEx$EXONEND,myEx$GENEID,max)
    
    olid1=olid2=NULL
    for (i in 1:(length(brStart)-1)){
      for (j in (i+1):length(brStart)){
        isParalog=FALSE
        myStatusij=(brStart[i]-brStart[j])/1e8*(brEnd[i]-brStart[j])
        myStatusji=(brStart[j]-brStart[i])/1e8*(brEnd[j]-brStart[i])
        if (myStatusij < 0 | myStatusji <0) isParalog=TRUE
        if (!isParalog){
          par1=c(as.character(mygene1[i]),geneParalog[which(geneParalog[,1]==as.character(mygene1[i])),2])
          par2=c(as.character(mygene1[j]),geneParalog[which(geneParalog[,1]==as.character(mygene1[j])),2])
          if (length(intersect(par1,par2))>0) isParalog=TRUE
        }
        if (isParalog){
          olid1=c(olid1,i)
          olid2=c(olid2,j)
        }
      }
    }
    if (length(olid1)>0){
      olGroup=seq(brStart)
      olGroup[olid2]=olid1
      olGroupU=rep(-1,length(olGroup))
      while(sum(olGroup!=olGroupU)>0){
        olGroupU=olGroup
        olGroup=olGroup[olGroup]
      }
      olGroupID=unique(olGroup)
      for (myol in olGroupID){#select only one
        keepID=which(olGroup==myol)
        if (length(keepID) >1){
          olFusionName=paste(mygene1[keepID],myge,sep="-")
          olFusion=dupFusion[match(olFusionName,dupFusion$name12),]
          selectID=which.max(olFusion$totalCount)
          #check if no fusion gene has higher total count
          if (sum(olFusion$totalCount[selectID] > olFusion$totalCount[-selectID])==0){
            #so we need another way to select the best candidate, we select the one closer geneDist
            selectID=which.min(olFusion$geneDist)
          }
          #update to the removal list
          keepID=which(!is.na(match(myFusion$name12,olFusionName)))
          myFusion$note5[keepID]=paste(myFusion$note5[keepID],paste(olFusionName,collapse = ","),sep=";")
        }
        
      }
    }
  }


##### checking read sequences

fastaDat2=fastaDat1=list();
fastaSplit2=fastaSplit1=list();
frfiles=list.files(inPath,paste(readStrands,"_fusionMappedReadsChunk_*",sep=""))
for (i in 1:length(frfiles)){
  #read fasta files of mapped reads
  ftag=rev(strsplit(strsplit(frfiles[i],"\\.")[[1]][1],"_")[[1]])[1]
  # con <- file(paste(inPath,"/",readStrands,"_fastaseq_",ftag,"_1.fa",sep=""), "r", blocking = FALSE)
  # mydata=readLines(con)
  # close(con)
  # fastaDat1[[frfiles[i]]]=mydata
  # 
  # con <- file(paste(inPath,"/",readStrands,"_fastaseq_",ftag,"_2.fa",sep=""), "r", blocking = FALSE)
  # mydata=readLines(con)
  # close(con)
  # fastaDat2[[frfiles[i]]]=mydata
  
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

myheader=as.character(myFusion$header)
rID.adj=match(myheader,fastaSplit1)
mySplit5N=fastaSplit1[rID.adj+1]
mySplit3N=fastaSplit2[rID.adj+1]
#if RF, convert 5N
if (readStrands=="RF" || readStrands=="UN" || readStrands=="RR") mySplit5N=sapply(mySplit5N,convertReverseComplement) else mySplit3N=sapply(mySplit3N,convertReverseComplement)


library(Biostrings)
fasta = readDNAStringSet(txFastaFile)
fasta_txnames=sapply(names(fasta), function(x) unlist(strsplit(x," "))[1])
#shrinkLen=5

tx5fa=rep("",length(mySplit5N))
frontTxList=as.character(myFusion$front_tx)
frontTxSet=unique(frontTxList)
cat("\n Total transcripts at 5 prime",length(frontTxSet))
for (i in 1:length(frontTxSet)){
  keepID=which(frontTxList==frontTxSet[i])
  txname=frontTxSet[i]
  matchID=match(txname,fasta_txnames)
  mytxFasta=fasta[matchID]
  mybrpos=myFusion$front_hitpos[keepID] -myFusion$front_querypos[keepID] +1
  tx5fa[keepID]=substring(mytxFasta,mybrpos,mybrpos+fragmentInfo$readlen-1)
}


tx3fa=rep("",length(mySplit3N))
backTxList=as.character(myFusion$back_tx)
backTxSet=unique(backTxList)
cat("\n Total transcripts at 3 prime",length(backTxSet))
for (i in 1:length(backTxSet)){
  keepID=which(backTxList==backTxSet[i])
  txname=backTxSet[i]
  matchID=match(txname,fasta_txnames)
  mytxFasta=fasta[matchID]
  mybrpos=myFusion$back_hitpos[keepID] - myFusion$back_querypos[keepID]+1
  tx3fa[keepID]=substring(mytxFasta,mybrpos,mybrpos+fragmentInfo$readlen-1)
}

#matchedNum.thres=fragmentInfo$readlen-fragmentInfo$kmer+1+10
matchedNum.thres=fragmentInfo$readlen*0.85
rmID=NULL
matchedNum=apply(cbind(mySplit5N,tx5fa),1,function(x){
  mync=min(nchar(x[1]),nchar(x[2]))
  return(sum(unlist(strsplit(x[1],""))[1:mync]==unlist(strsplit(x[2],""))[1:mync]))
})
rmID=c(rmID,which(matchedNum>matchedNum.thres))
myFusion$mn5N5=matchedNum

matchedNum=apply(cbind(mySplit3N,tx5fa),1,function(x){
  mync=min(nchar(x[1]),nchar(x[2]))
  return(sum(unlist(strsplit(x[1],""))[1:mync]==unlist(strsplit(x[2],""))[1:mync]))
})
rmID=c(rmID,which(matchedNum>matchedNum.thres))
myFusion$mn3N5=matchedNum

matchedNum=apply(cbind(mySplit5N,tx3fa),1,function(x){
  mync=min(nchar(x[1]),nchar(x[2]))
  return(sum(unlist(strsplit(x[1],""))[(nchar(x[1])-mync+1):nchar(x[1])]==unlist(strsplit(x[2],""))[(nchar(x[2])-mync+1):nchar(x[2])]))
})
rmID=c(rmID,which(matchedNum>matchedNum.thres))
myFusion$mn5N3=matchedNum

matchedNum=apply(cbind(mySplit3N,tx3fa),1,function(x){
  mync=min(nchar(x[1]),nchar(x[2]))
  return(sum(unlist(strsplit(x[1],""))[(nchar(x[1])-mync+1):nchar(x[1])]==unlist(strsplit(x[2],""))[(nchar(x[2])-mync+1):nchar(x[2])]))
})
rmID=c(rmID,which(matchedNum>matchedNum.thres))
myFusion$mn3N3=matchedNum

rmID=unique(rmID)
if (length(rmID)>0) myFusion=myFusion[-rmID,]


#####check the direction of read pairs:  
#FR: correct direction is 1-0 and 2-1
#RF or UN: correct direction is 2-0 and 1-1
if (readStrands=="RF" || readStrands=="RR")  rmID=c(which(myFusion$readType==1 & myFusion$direction==0),which(myFusion$readType==2 & myFusion$direction==1))

if (readStrands=="FR" || readStrands=="FF")  rmID=c(which(myFusion$readType==1 & myFusion$direction==1),which(myFusion$readType==2 & myFusion$direction==0))

if (readStrands=="UN") rmID=NULL #can not apply

#remove the fusion candidates if there are at least 5% incorrect direction read pairs
if (length(rmID)>0){
  res=table(myFusion$name12[rmID])
  myFusion$fDirCount=0
  myFusion$fDirCount[which(!is.na(match(as.character(myFusion$name12),names(res))))]=res[na.omit(match(as.character(myFusion$name12),names(res)))]
  myFusion$fDirProp=myFusion$fDirCount/myFusion$supportCount
  myFusion=myFusion[myFusion$fDirProp < 0.05,]
}



#### filter by shared break points when the partner genes are not paralogs
myFusion2=myFusion
res=table(myFusion2$name12)
myFusion2$supportCount3=res[match(myFusion2$name12,names(res))]
myFusion2=myFusion2[myFusion2$supportCount3 > 1,] #consider at least 2 supporting split reads
myDup=duplicated(myFusion2$name12)
myFusionNoDup=myFusion2[!myDup,]
dupge1=table(myFusionNoDup$gene1)
dupge2=table(myFusionNoDup$gene2)
myFusion2=cbind(myFusion2,as.integer(dupge1[match(as.character(myFusion2$gene1),names(dupge1))]))
myFusion2=cbind(myFusion2,as.integer(dupge2[match(as.character(myFusion2$gene2),names(dupge2))]))
colnames(myFusion2)[c(ncol(myFusion2)-1,ncol(myFusion2))]=c("dupge1f2","dupge2f2")

frontDupge=unique(as.character(myFusion2$gene1[which(myFusion2$dupge1f2>1)]))
rmFusion=NULL;
if (length(frontDupge)>0)
  for (myge in frontDupge){
    dupFusion=myFusion2[which(as.character(myFusion2$gene1)==myge),]
    dupFusion=dupFusion[!duplicated(dupFusion$name12),]
    #group by brchposEx5
    brPos=dupFusion$brchposEx5
    names(brPos)=seq(brPos)
    brPos=sort(brPos, decreasing = FALSE)
    myGroup=seq(brPos)
    myDiff=diff(brPos)
    myGroup[which(myDiff<=shrinkLen)+1]=which(myDiff<=shrinkLen)
    myGroupU=rep(-1,length(myGroup))
    while(sum(myGroup!=myGroupU)>0){
      myGroupU=myGroup
      myGroup=myGroup[myGroup]
    }
    myGroupID=unique(myGroup)
    for (myGrp in myGroupID){
      keepID=which(myGroup==myGrp)
      if (length(keepID) >1){
        grpFusion=dupFusion[as.integer(names(brPos[keepID])),]
        if(sum(grpFusion$note3=="" & grpFusion$note5=="") >0) rmFusion=c(rmFusion,unique(grpFusion$name12))
      }
    }
    
  }

backDupge=unique(as.character(myFusion2$gene2[which(myFusion2$dupge2f2>1)]))
if (length(backDupge)>0)
  for (myge in backDupge){
    dupFusion=myFusion2[which(as.character(myFusion2$gene2)==myge),]
    dupFusion=dupFusion[!duplicated(dupFusion$name12),]
    #group by brchposEx3
    brPos=dupFusion$brchposEx3
    names(brPos)=seq(brPos)
    brPos=sort(brPos, decreasing = FALSE)
    myGroup=seq(brPos)
    myDiff=diff(brPos)
    myGroup[which(myDiff<=shrinkLen)+1]=which(myDiff<=shrinkLen)
    myGroupU=rep(-1,length(myGroup))
    while(sum(myGroup!=myGroupU)>0){
      myGroupU=myGroup
      myGroup=myGroup[myGroup]
    }
    myGroupID=unique(myGroup)
    for (myGrp in myGroupID){
      keepID=which(myGroup==myGrp)
      if (length(keepID) >1){
        grpFusion=dupFusion[as.integer(names(brPos[keepID])),]
        if(sum(grpFusion$note3=="" & grpFusion$note5=="") >0) rmFusion=c(rmFusion,unique(grpFusion$name12))
      }
    }
  }




if (length(rmFusion)){
  rmID=which(!is.na(match(as.character(myFusion$name12),rmFusion)))
  myFusion=myFusion[-rmID,]
}


# filter by junction distances
rmID=which(as.character(myFusion$chrom1)==as.character(myFusion$chrom2) & myFusion$junctDist<=FuSeq.params$minJunctionDist)
if (length(rmID)>0) myFusion=myFusion[-rmID,]

res=table(myFusion$tx12)
myFusion$tx12Count3=res[match(myFusion$tx12,names(res))]
myFusion$mEstCount=myFusion$tx12Count3/(fragmentInfo$readlen-2*fragmentInfo$kmer-1) * (2*fragmentInfo$kmer+1)
myFusion$mEstProp=myFusion$mappedCount/myFusion$mEstCount



############

#check the consistency between mapped reads and split reads
#myFusionMapped=myFusionTmp2[myFusionTmp2$mappedCrtCount>=1,]
myFusion$fusionName=myFusion$name12
myFusionMapped=myFusion
myDup=duplicated(myFusionMapped$name12)
myFusionMapped=myFusionMapped[!myDup,]


if (readStrands=="RF" || readStrands=="UN" || readStrands=="RR"){
  #switch name21 - name12 for running function detectJunctionBreaks() 
  myFusionMapped$name12_raw=myFusionMapped$name12
  myFusionMapped$name21_raw=myFusionMapped$name21
  myFusionMapped$name12=myFusionMapped$name21_raw
  myFusionMapped$name21=myFusionMapped$name12_raw
  
  myFusionMapped$gene1_raw=myFusionMapped$gene1
  myFusionMapped$gene2_raw=myFusionMapped$gene2
  myFusionMapped$gene1=myFusionMapped$gene2_raw
  myFusionMapped$gene2=myFusionMapped$gene1_raw
}

matchID=match(as.character(myFusionMapped$name12), names(FuSeq.MR$feqInfo$feqFgeMap))
myFusionMapped=myFusionMapped[which(!is.na(matchID)),]

if (nrow(myFusionMapped) > 0){
  junctBr=detectJunctionBreaks(myFusionMapped,inPath, FuSeq.MR$feqInfo$feq,FuSeq.MR$feqInfo$feqFgeMap, anntxdb, readStrands=FuSeq.params$readStrands)
  myFusionMapped=junctBr$myFusionFinal
  
  ok5Pos=ok3Pos=NULL
  for (i in 1:nrow(myFusionMapped)){
    res=junctBr$junctInfo[[myFusionMapped$fusionName[i]]]
    ck5Pos=unlist(res$readL.chrPos)
    if (myFusionMapped$strand1[i]=="+") ck5Pos= sum(ck5Pos+unlist(res$readL.seqLen)-1-shrinkLen<= myFusionMapped$brchposEx5[i])/length(ck5Pos) else ck5Pos= sum(ck5Pos-unlist(res$readL.seqLen)+1+shrinkLen >= myFusionMapped$brchposEx5[i])/length(ck5Pos)
    ok5Pos=c(ok5Pos,ck5Pos)
    
    ck3Pos=unlist(res$readR.chrPos)
    if (myFusionMapped$strand2[i]=="+") ck3Pos= sum(ck3Pos+shrinkLen>= myFusionMapped$brchposEx3[i])/length(ck3Pos) else ck3Pos= sum(ck3Pos-shrinkLen<= myFusionMapped$brchposEx3[i])/length(ck3Pos)
    ok3Pos=c(ok3Pos,ck3Pos)
  }
  myFusionMapped$ok5Pos=ok5Pos
  myFusionMapped$ok3Pos=ok3Pos
  
  rmFusion=c(myFusionMapped$name12[myFusionMapped$ok5Pos==0],myFusionMapped$name12[myFusionMapped$ok3Pos==0])

  #txlen test
  tx3LenTest=tx5LenTest=NULL;
  for (i in 1:nrow(myFusionMapped)){
    res=testFtxlen(mu=fragmentInfo$fragLengthMean, sig=fragmentInfo$fragLengthSd, r=fragmentInfo$readlen,ftxlen=myFusionMapped$flen5[i], kmerlen=fragmentInfo$kmer,fragDist=fragDist,M=10000)
    tx5LenTest=c(tx5LenTest,sum(res$x<=myFusionMapped$tx5LenMean[i])/length(res$x))
    res=testFtxlen(mu=fragmentInfo$fragLengthMean, sig=fragmentInfo$fragLengthSd, r=fragmentInfo$readlen,ftxlen=myFusionMapped$flen3[i], kmerlen=fragmentInfo$kmer,fragDist=fragDist,M=10000)
    tx3LenTest=c(tx3LenTest,sum(res$x<=myFusionMapped$tx3LenMean[i])/length(res$x))
  }
  myFusionMapped$tx5LenTest=tx5LenTest
  myFusionMapped$tx3LenTest=tx3LenTest
  rmID=c(which(myFusionMapped$tx5LenSd<=1),which(myFusionMapped$tx3LenSd<=1),which(myFusionMapped$tx5LenTest<=0.05),which(myFusionMapped$tx3LenTest<=0.05))

  #estimate split reads from mapped reads
  myFusionMapped$srEst3=myFusionMapped$fragmentNum/myFusionMapped$flen3 * (fragmentInfo$readlen-2*fragmentInfo$kmer)
  myFusionMapped$srEst5=myFusionMapped$fragmentNum/myFusionMapped$flen5 * (fragmentInfo$readlen-2*fragmentInfo$kmer)
  myFusionMapped$srEstCount=apply(cbind(myFusionMapped$srEst3,myFusionMapped$srEst5),1,max)

  #junctInfo=junctBr$junctInfo
  exonInfo=select(anntxdb, keys=unique(c(as.character(myFusionMapped$front_tx),as.character(myFusionMapped$back_tx))), columns=c("EXONSTART","EXONEND"), keytype = "TXNAME")
  myFusionMapped$tx2End=myFusionMapped$tx2Start=myFusionMapped$tx1End=myFusionMapped$tx1Start=NULL
  myFusionMapped$tx2exAN=myFusionMapped$tx1exAN=NULL
  for (i in 1:nrow(myFusionMapped)){
    txname=as.character(myFusionMapped$front_tx[i])
    res=exonInfo[exonInfo$TXNAME==txname,]
    myFusionMapped$tx1Start[i]=min(res$EXONSTART)
    myFusionMapped$tx1End[i]=max(res$EXONEND)
    if (myFusionMapped$strand1[i]=="+") myFusionMapped$tx1exAN[i]=sum(res$EXONSTART<= myFusionMapped$brpos5.start[i]) else myFusionMapped$tx1exAN[i]=sum(res$EXONEND>= myFusionMapped$brpos5.start[i])
    
    txname=as.character(myFusionMapped$back_tx[i])
    res=exonInfo[exonInfo$TXNAME==txname,]
    myFusionMapped$tx2Start[i]=min(res$EXONSTART)
    myFusionMapped$tx2End[i]=max(res$EXONEND)
    if (myFusionMapped$strand2[i]=="+") myFusionMapped$tx2exAN[i]=sum(res$EXONEND>= myFusionMapped$brpos3.start[i]) else  myFusionMapped$tx2exAN[i]=sum(res$EXONSTART<= myFusionMapped$brpos3.start[i])
    
  }
  
  
  exonInfo=select(anntxdb, keys=unique(c(as.character(myFusionMapped$gene1),as.character(myFusionMapped$gene2))), columns=c("GENEID","TXNAME","EXONID","EXONSTART","EXONEND","TXSTRAND"), keytype = "GENEID")
  myFusionMapped$gechr.start1=myFusionMapped$gechr.end1=myFusionMapped$exid.start1=myFusionMapped$exid.end1=NULL
  myFusionMapped$gechr.start2=myFusionMapped$gechr.end2=myFusionMapped$exid.start2myFusionMapped$exid.end2=NULL
  for (i in 1:nrow(myFusionMapped)){
    myEx=exonInfo[exonInfo$GENEID==myFusionMapped$gene1[i],]
    myFusionMapped$gechr.start1[i]=min(myEx$EXONSTART)
    myFusionMapped$gechr.end1[i]=max(myEx$EXONEND)
    myFusionMapped$exid.start1[i]=myEx$EXONID[which.min(myEx$EXONSTART)]
    myFusionMapped$exid.end1[i]=myEx$EXONID[which.max(myEx$EXONEND)]
  
    myEx=exonInfo[exonInfo$GENEID==myFusionMapped$gene2[i],]
    myFusionMapped$gechr.start2[i]=min(myEx$EXONSTART)
    myFusionMapped$gechr.end2[i]=max(myEx$EXONEND)
    myFusionMapped$exid.start2[i]=myEx$EXONID[which.min(myEx$EXONSTART)]
    myFusionMapped$exid.end2[i]=myEx$EXONID[which.max(myEx$EXONEND)]
    
  }
  
  #extract junctInfo
  junctInfo=junctBr$junctInfo
  
  exonInfo=select(anntxdb, keys=unique(c(as.character(myFusionMapped$gene1),as.character(myFusionMapped$gene2))), columns=c("GENEID","TXNAME","EXONID","EXONSTART","EXONEND","TXSTRAND"), keytype = "GENEID")
  for (i in 1:nrow(myFusionMapped)){
    res=junctInfo[[i]]

    mappedExonL=exonInfo[match(sort(unique(unlist(res$readL.ExonID))),exonInfo$EXONID),]
    geneExonL=exonInfo[exonInfo$GENEID==mappedExonL$GENEID[1],]
    
    mappedExonR=exonInfo[match(sort(unique(unlist(res$readR.ExonID))),exonInfo$EXONID),]
    geneExonR=exonInfo[exonInfo$GENEID==mappedExonR$GENEID[1],]
  
    junctInfo[[i]]$mappedExonL=mappedExonL
    junctInfo[[i]]$geneExonL=geneExonL
    
    junctInfo[[i]]$mappedExonR=mappedExonR
    junctInfo[[i]]$geneExonR=geneExonR
  }
  
 
  exonInfo=select(anntxdb, keys=unique(c(as.character(myFusionMapped$gene1),as.character(myFusionMapped$gene2))), columns=c("GENEID","TXNAME","EXONID","EXONSTART","EXONEND","TXSTRAND"), keytype = "GENEID")
  myFusionMapped$ANflen3=myFusionMapped$ANflen5=NULL
  for (i in 1:nrow(myFusionMapped)){
    res=junctInfo[[i]]
    #summary(unlist(res$new.readL.GenePos))
    
    mappedExonL=exonInfo[match(sort(unique(unlist(res$readL.ExonID))),exonInfo$EXONID),]
    geneExonL=exonInfo[exonInfo$GENEID==mappedExonL$GENEID[1],]
    
    mappedExonR=exonInfo[match(sort(unique(unlist(res$readR.ExonID))),exonInfo$EXONID),]
    geneExonR=exonInfo[exonInfo$GENEID==mappedExonR$GENEID[1],]
    
    ###get annotated flen
    if (mappedExonL$TXSTRAND[1]=="+") ANExonL=geneExonL[geneExonL$EXONSTART <= min(mappedExonL$EXONSTART),] else ANExonL=geneExonL[geneExonL$EXONSTART >= min(mappedExonL$EXONSTART),]
      myFusionMapped$ANflen5[i]=getGeneLen(geneExonMat=ANExonL)
      if (mappedExonR$TXSTRAND[1]=="+") ANExonR=geneExonR[geneExonR$EXONSTART >= min(mappedExonR$EXONSTART),] else ANExonR=geneExonR[geneExonR$EXONSTART <= min(mappedExonR$EXONSTART),]
      myFusionMapped$ANflen3[i]=getGeneLen(geneExonMat=ANExonR)
      

  #Estimate split reads from mapped reads
  myFusion$srEstCount=0
  myFusion$srEstCount[which(!is.na(match(as.character(myFusion$fusionName),as.character(myFusionMapped$fusionName))))]=myFusionMapped$srEstCount[na.omit(match(as.character(myFusion$fusionName),as.character(myFusionMapped$fusionName)))]
  # filter by the estimated number of split reads from mapped reads
  rmID=which(myFusion$mappedCount>0 & myFusion$supportCount < trunc(myFusion$srEstCount)-1)
  if (length(rmID)>0) myFusion=myFusion[-rmID,]
  
  #filter by standard deviation of flen of the mapped reads
  rmID=c(which(myFusionMapped$tx5LenSd<=1 & myFusionMapped$fragmentNum>=3),which(myFusionMapped$tx3LenSd<=1 & myFusionMapped$fragmentNum>=3))
  rmFusion=unique(c(myFusionMapped$fusionName[rmID]))
  if (length(rmFusion)>0){
    myFusion=myFusion[which(is.na(match(myFusion$fusionName,rmFusion))),]
  }
  }
  
  ### consistency between mapped and split - we compare by the end points to be more robust
  myFusionMapped$brsign5=ifelse(myFusionMapped$brchposEx5-myFusionMapped$brpos5.end>0,"+","-")
  myFusionMapped$brSign3=ifelse(myFusionMapped$brchposEx3-myFusionMapped$brpos3.end>0,"-","+")
  # This filter depends very much on the correction of the mapped reads. In most of cases, they are ok. However, if only 1 wrong mapped read pair might change the results. We use this filter as a soft filter after this step.
  rmID=c(which(myFusionMapped$strand1!=myFusionMapped$brsign5),which(myFusionMapped$strand2!=myFusionMapped$brsign3))
  if (length(rmID)>0){
    rmFusion=myFusionMapped$fusionName[rmID]
    myFusion=myFusion[which(is.na(match(myFusion$fusionName,rmFusion))),]
  }
  
  junctBr.refine=refineJunctionBreak(junctBr, anntxdb, fragmentInfo, readStrands, fragDist)
  
} else {
  myFusionMapped=NULL
  junctBr.refine=NULL
  junctInfo=NULL
} #end of if nrow(myFusionMapped)>0

### proportion of remaining reads after filtering
res=table(myFusion$name12)
myFusion$supportCount3=res[match(myFusion$name12,names(res))]
myFusion$survProp=myFusion$supportCount3/myFusion$supportCount

### check for fusion genes with supporting reads >=5
myFusionFinal=myFusion
myFusionFinal=myFusionFinal[!duplicated(myFusionFinal$name12),]
dim(myFusionFinal)
keepID=which(myFusionFinal$supportCount3>=5)
rmID=NULL
if (length(keepID)>0)
for (i in 1:length(keepID)){
  myID=which(as.character(myFusion$name12)==as.character(myFusionFinal$name12[keepID[i]]))
  res=myFusion[myID,]
  
  if(sd(res$front_hitpos)<1)  rmID=c(rmID,myID)
  if(sd(res$back_hitpos)<1)  rmID=c(rmID,myID)
}

if(length(rmID)>0) myFusion=myFusion[-rmID,]



### get final results
myFusionFinal=myFusion
myFusionFinal=myFusionFinal[!duplicated(myFusionFinal$name12),]
dim(myFusionFinal)

# Filter by minSR
myFusionFinal=myFusionFinal[myFusionFinal$supportCount>=FuSeq.params$minSR,]
# Check if fusion breaks locate at the ending exon of gene
rmID.SR=checkEndExon(myFusionFinal=myFusionFinal, junctBr=NULL, anntxdb=anntxdb, readStrands=FuSeq.params$readStrands,shrinkLen=5, type="SR")

if (length(rmID.SR)>0) myFusionFinal=myFusionFinal[-rmID.SR,]

#assign a score by totalCount
myFusionFinal$score=myFusionFinal$totalCount

myFusionFinal$supportRead=myFusionFinal$supportCount+myFusionFinal$mappedCount

return(list(junctInfo=junctInfo,myFusionFinal=myFusionFinal, myFusion=myFusion, myFusionMapped=myFusionMapped,junctBr.refine=junctBr.refine))


}




