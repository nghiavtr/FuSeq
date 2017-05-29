############################################################
#####post process mapped reads
############################################################

postProcessMappedRead <-function(inPath, anntxdb, FuSeq.SR, FuSeq.MR, FuSeq.params){
  cat("\n Post processing mapped reads (MR)...")
  
  myFusionFinal=FuSeq.MR$myFusionFinal
  fragmentInfo=FuSeq.SR$fragmentInfo
  fragDist=FuSeq.SR$fragDist
  
  ##### revisit the information of junction breaks
  junctBr=refineJunctionBreak(FuSeq.MR$junctBr, anntxdb, fragmentInfo, FuSeq.params$readStrands, fragDist=fragDist)
  # update to myFusionFinal
  
  myFusionFinal$adjFragmentNum=myFusionFinal$outlierNum=NA
  myFusionFinal$ftxLenProp=myFusionFinal$ftxMedianLen=myFusionFinal$ftxMeanLen=NA
  myFusionFinal$ftx5LenSd=myFusionFinal$ftx5LenMean=myFusionFinal$ftx3LenSd=myFusionFinal$ftx3LenMean=NA
  
  matchID=match(as.character(myFusionFinal$fusionName),as.character(junctBr$myFusionFinal$fusionName))
  
  myFusionFinal$adjFragmentNum=junctBr$myFusionFinal$adjFragmentNum[matchID]
  myFusionFinal$outlierNum= junctBr$myFusionFinal$outlierNum[matchID]
  myFusionFinal$ftxLenProp=junctBr$myFusionFinal$ftxLenProp[matchID]
  myFusionFinal$ftxMedianLen=junctBr$myFusionFinal$ftxMedianLen[matchID]
  myFusionFinal$ftxMeanLen=junctBr$myFusionFinal$ftxMeanLen[matchID]
  myFusionFinal$ftx5LenMean=junctBr$myFusionFinal$ftx5LenMean[matchID]
  myFusionFinal$ftx5LenSd=junctBr$myFusionFinal$ftx5LenSd[matchID]
  myFusionFinal$ftx3LenMean=junctBr$myFusionFinal$ftx3LenMean[matchID]
  myFusionFinal$ftx3LenSd=junctBr$myFusionFinal$ftx3LenSd[matchID]
  
  keepID=which(!is.na(junctBr$myFusionFinal$crt.tx5LenMean[matchID]))
  if (length(keepID)>0){
    myFusionFinal$tx5LenMean[keepID]=junctBr$myFusionFinal$crt.tx5LenMean[matchID][keepID]
    myFusionFinal$tx5LenSd[keepID]=junctBr$myFusionFinal$crt.tx5LenSd[matchID][keepID]
    myFusionFinal$tx3LenMean[keepID]=junctBr$myFusionFinal$crt.tx3LenMean[matchID][keepID]
    myFusionFinal$tx3LenSd[keepID]=junctBr$myFusionFinal$crt.tx3LenSd[matchID][keepID]
  }
  
  ##### quick filter here
  keepID=which(myFusionFinal$ftx5LenMean > 0)
  myFusionFinal=myFusionFinal[keepID,]
  keepID=which(myFusionFinal$ftx3LenMean > 0)
  myFusionFinal=myFusionFinal[keepID,]
  
  keepID=which(myFusionFinal$ftxMeanLen < 1000 & myFusionFinal$ftxMeanLen > fragmentInfo$readlen )
  myFusionFinal=myFusionFinal[keepID,]
  #get empirical probabilities
  myFusionFinal$ftxLenEmpProp=sapply(myFusionFinal$ftxMedianLen,function(x) sum(fragDist[fragDist[,1]>x,2])/sum(fragDist[,2]))
  #filter by 1e-04
  myFusionFinal=myFusionFinal[myFusionFinal$ftxLenEmpProp>1e-04,]
  
  
  ########################## filtering
  #Remove false positives
  matchID=match(myFusionFinal$fusionName,FuSeq.SR$myFusionFP$name12)
  myFusionFinal=myFusionFinal[which(is.na(matchID)),]
  #txlen test
  rmID=c(which(myFusionFinal$tx5LenSd<=1 & myFusionFinal$fragmentNum>=3),which(myFusionFinal$tx3LenSd<=1 & myFusionFinal$fragmentNum>=3))
  rmFusion=unique(c(myFusionFinal$fusionName[rmID]))
  if (length(rmFusion)>0){
    myFusionFinal=myFusionFinal[which(is.na(match(myFusionFinal$fusionName,rmFusion))),]
  }
  
  
  tx3LenTest=tx5LenTest=NULL;
  for (i in 1:nrow(myFusionFinal)){
    res=testFtxlen(mu=fragmentInfo$fragLengthMean, sig=fragmentInfo$fragLengthSd, r=fragmentInfo$readlen,ftxlen=myFusionFinal$flen5[i], kmerlen=fragmentInfo$kmer,fragDist=fragDist,M=10000)
    tx5LenTest=c(tx5LenTest,sum(res$x<=myFusionFinal$tx5LenMean[i])/length(res$x))
    res=testFtxlen(mu=fragmentInfo$fragLengthMean, sig=fragmentInfo$fragLengthSd, r=fragmentInfo$readlen,ftxlen=myFusionFinal$flen3[i], kmerlen=fragmentInfo$kmer,fragDist=fragDist,M=10000)
    tx3LenTest=c(tx3LenTest,sum(res$x<=myFusionFinal$tx3LenMean[i])/length(res$x))
  }
  myFusionFinal$tx5LenTest=tx5LenTest
  myFusionFinal$tx3LenTest=tx3LenTest
  myFusionFinal=myFusionFinal[myFusionFinal$tx5LenTest >= 0.10,]
  myFusionFinal=myFusionFinal[myFusionFinal$tx3LenTest >= 0.10,]
  
  # ##### create a score
  myFusionFinal$score=myFusionFinal$correctedCount
  
  #get split reads  
  myFusionFinal$SR_supportCount=0
  matchID=match(myFusionFinal$fusionName,FuSeq.SR$fusionGene$name12)
  myFusionFinal$SR_supportCount[which(!is.na(matchID))]=FuSeq.SR$fusionGene$supportCount[na.omit(matchID)]
  
  myFusionFinal$SR_correctedCount=0
  matchID=match(myFusionFinal$fusionName,FuSeq.SR$fusionGene$name12)
  myFusionFinal$SR_correctedCount[which(!is.na(matchID))]=FuSeq.SR$fusionGene$adjtx12Count[na.omit(matchID)]
  myFusionFinal$score=myFusionFinal$score+myFusionFinal$SR_correctedCount
  
  
  myFusionFinal$srEst3=myFusionFinal$supportCount/myFusionFinal$flen3 * (fragmentInfo$readlen-2*fragmentInfo$kmer)
  myFusionFinal$srEst5=myFusionFinal$supportCount/myFusionFinal$flen5 * (fragmentInfo$readlen-2*fragmentInfo$kmer)
  #smRatio indicates the ratio between mapped reads and split reads, we expect 1 split reads randomly happen from 1000 mapped reads
  smRatio=0.001
  rmID=unique(c(which(myFusionFinal$SR_supportCount > smRatio*myFusionFinal$supportCount & myFusionFinal$SR_supportCount<trunc(myFusionFinal$srEst5)),which(myFusionFinal$SR_supportCount>smRatio*myFusionFinal$supportCount & myFusionFinal$SR_supportCount<trunc(myFusionFinal$srEst3))))
  if (length(rmID)>0) myFusionFinal=myFusionFinal[-rmID,]
  
  #check ending exon in MR
  rmID.MR=checkEndExon(myFusionFinal=myFusionFinal, junctBr=junctBr, anntxdb=anntxdb, readStrands=FuSeq.params$readStrands,shrinkLen=5, type="MR")
  #cat("\n Number of ending exon in MR: ",length(rmID.MR))
  if (length(rmID.MR)>0) myFusionFinal=myFusionFinal[-rmID.MR,]
  
  myFusionFinal$supportRead=myFusionFinal$supportCount+myFusionFinal$SR_supportCount
  
  res=list(myFusionFinal=myFusionFinal, junctBr.refine=junctBr)
  return(res)
}



