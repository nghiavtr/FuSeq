############################################################
### 08 November 2018
# - Ignore if some biological databases, e.g.,ribSubunitDb,mitoTransDb,ribonuproDb do not exist
# - Add information of exon boundary that matches to breaking points
##### Combine two fusion gene candidate lists
integrateFusion <-function(myFusionFinal.MR, myFusionFinal.SR, FuSeq.params, fragmentInfo=NULL, paralog.fc.thres=NULL){

  minScore=FuSeq.params$minScore
  #####last filter for MR
  #filter by score
  myFusionFinal.MR=myFusionFinal.MR[myFusionFinal.MR$score>=minScore,]
  #junctiondistance
  rmID=which(as.character(myFusionFinal.MR$chrom1)==as.character(myFusionFinal.MR$chrom2) & myFusionFinal.MR$junctDist<=FuSeq.params$minJunctionDist)
  if (length(rmID)>0) myFusionFinal.MR=myFusionFinal.MR[-rmID,]
  #genedistance
  rmID=which(as.character(myFusionFinal.MR$chrom1)==as.character(myFusionFinal.MR$chrom2) & myFusionFinal.MR$geneDist<=FuSeq.params$minGeneDist)
  if (length(rmID)>0) myFusionFinal.MR=myFusionFinal.MR[-rmID,]
  
  
  ##### last filter for SR
  #junctionDistance
  rmID=which(as.character(myFusionFinal.SR$chrom1)==as.character(myFusionFinal.SR$chrom2) & myFusionFinal.SR$junctDist<=FuSeq.params$minJunctionDist)
  if (length(rmID)>0) myFusionFinal.SR=myFusionFinal.SR[-rmID,]
  # #geneDistance
  rmID=which(as.character(myFusionFinal.SR$chrom1)==as.character(myFusionFinal.SR$chrom2) & myFusionFinal.SR$geneDist<=FuSeq.params$minGeneDist)
  if (length(rmID)>0) myFusionFinal.SR=myFusionFinal.SR[-rmID,]
  
  #filter by score
  myFusionFinal.SR=myFusionFinal.SR[myFusionFinal.SR$totalCount>=minScore,]
  #keep only one presentative for each fusion gene - no, should allow multiple break points
  #myFusionFinal.SR=myFusionFinal.SR[!duplicated(myFusionFinal.SR$name12),]
  
  #dup check here
  myDup=duplicated(myFusionFinal.SR$name12)
  myFusionNoDup=myFusionFinal.SR[!myDup,]
  dupge1=table(myFusionNoDup$gene1)
  dupge2=table(myFusionNoDup$gene2)
  myFusionFinal.SR=cbind(myFusionFinal.SR,as.integer(dupge1[match(as.character(myFusionFinal.SR$gene1),names(dupge1))]))
  myFusionFinal.SR=cbind(myFusionFinal.SR,as.integer(dupge2[match(as.character(myFusionFinal.SR$gene2),names(dupge2))]))
  colnames(myFusionFinal.SR)[c(ncol(myFusionFinal.SR)-1,ncol(myFusionFinal.SR))]=c("dupge1f2","dupge2f2")
  rmID=which(myFusionFinal.SR$dupge1f2>1 & myFusionFinal.SR$dupge2f2>1)
  if (length(rmID)>0) myFusionFinal.SR=myFusionFinal.SR[-rmID,]
  
  
  #deal with paralogs
  if (!is.null(paralog.fc.thres)){
    rmFusionName=NULL
    res5=table(myFusionFinal.SR$note5)
    res5=res5[-which(names(res5)=="")]
    if (length(res5)>0){
      for (i in 1:length(res5)){
        mySR=myFusionFinal.SR[myFusionFinal.SR$note5==names(res5)[i],]
        mySR=mySR[order(mySR$totalCount,decreasing = TRUE),]
        rmID=which(mySR$totalCount*paralog.fc.thres<mySR$totalCount[1])
        if (length(rmID)>0){
          rmFusionName=c(rmFusionName,mySR$fusionName[rmID])
        }
      }
    }
    res3=table(myFusionFinal.SR$note3)
    res3=res3[-which(names(res3)=="")]
    if (length(res3)>0){
      for (i in 1:length(res3)){
        mySR=myFusionFinal.SR[myFusionFinal.SR$note3==names(res3)[i],]
        mySR=mySR[order(mySR$totalCount,decreasing = TRUE),]
        rmID=which(mySR$totalCount*paralog.fc.thres<mySR$totalCount[1])
        if (length(rmID)>0){
          rmFusionName=c(rmFusionName,mySR$fusionName[rmID])
        }
      }
    }
    
    if (length(rmFusionName)>0){
      rmFusionName=unique(rmFusionName)
      rmID=match(rmFusionName,myFusionFinal.SR$fusionName)
      myFusionFinal.SR=myFusionFinal.SR[-rmID,]
    }
  }
  #if the fusion shares breaking points in both sides
  rmID=which(myFusionFinal.SR$note3!="" & myFusionFinal.SR$note5!="")
  if (length(rmID)>0) myFusionFinal.SR=myFusionFinal.SR[-rmID,]

  # create common features
  myFusionFinal.SR$chrom5p=myFusionFinal.SR$chrom1
  myFusionFinal.SR$brpos5.start=myFusionFinal.SR$brchposEx5
  myFusionFinal.SR$brpos5.end=myFusionFinal.SR$brchposEx5+1
  myFusionFinal.SR$strand5p=myFusionFinal.SR$strand1
  myFusionFinal.SR$chrom3p=myFusionFinal.SR$chrom2
  myFusionFinal.SR$brpos3.start=myFusionFinal.SR$brchposEx3
  myFusionFinal.SR$brpos3.end=myFusionFinal.SR$brchposEx3+1
  myFusionFinal.SR$strand3p=myFusionFinal.SR$strand2

  #make columns first
  myFusionFinal.SR$exonID5p= myFusionFinal.SR$brpos5.start
  myFusionFinal.SR$exonBound5p= myFusionFinal.SR$brpos5.start
  myFusionFinal.SR$exonID3p= myFusionFinal.SR$brpos5.start
  myFusionFinal.SR$exonBound3p= myFusionFinal.SR$brpos5.start


####### find exon boundary
### re-adjust breaking points in split reads: due to the similarity between the end-sequence of 5p and the start sequence of 3p, the breaking points can be shifted a abit, specially for fusions with two genes from different chromosomes.
  # - 5p gene: for plus strand position is generally matched to exon-end; and exon-start for minus strand
  # - 3p gene: for plus strand position is generally matched to exon-start; and exon-end for minus strand
  if (nrow(myFusionFinal.SR)>0){
    myFusionFinal.SR$exonID5p=-1
    myFusionFinal.SR$exonBound5p=-1
    myFusionFinal.SR$exonID3p=-1
    myFusionFinal.SR$exonBound3p=-1
    if (nrow(myFusionFinal.SR)>0){
      exonInfo=select(anntxdb, keys=unique(c(as.character(myFusionFinal.SR$front_tx),as.character(myFusionFinal.SR$back_tx))), columns=c("TXNAME","TXSTRAND","EXONID","EXONSTART","EXONEND"), keytype = "TXNAME")
      exonBound5p=exonBound3p=NULL
      exonID5p=exonID3p=NULL
      for (i in 1:nrow(myFusionFinal.SR)){
        x=myFusionFinal.SR[i,]  
        br5=x$brchposEx5
        tx5=as.character(x$front_tx)
        e5=exonInfo[exonInfo$TXNAME == tx5,]
        if (e5$TXSTRAND[1]=="+"){
          absdif=abs(br5-e5$EXONEND)
          k=which.min(absdif)
          exonID5p=c(exonID5p,e5$EXONID[k])
          exonBound5p=c(exonBound5p,e5$EXONEND[k])
        }
        if (e5$TXSTRAND[1]=="-"){
          absdif=abs(br5-e5$EXONSTART)
          k=which.min(absdif)
          exonID5p=c(exonID5p,e5$EXONID[k])
          exonBound5p=c(exonBound5p,e5$EXONSTART[k])    
        }
        br3=x$brchposEx3
        tx3=as.character(x$back_tx)
        e3=exonInfo[exonInfo$TXNAME == tx3,]
        if (e3$TXSTRAND[1]=="+"){
          absdif=abs(br3-e3$EXONSTART)
          k=which.min(absdif)
          exonID3p=c(exonID3p,e3$EXONID[k])
          exonBound3p=c(exonBound3p,e3$EXONSTART[k])
        }
        if (e3$TXSTRAND[1]=="-"){
          absdif=abs(br3-e3$EXONEND)
          k=which.min(absdif)
          exonID3p=c(exonID3p,e3$EXONID[k])
          exonBound3p=c(exonBound3p,e3$EXONEND[k])
        }
      }
      myFusionFinal.SR$exonID5p=exonID5p
      myFusionFinal.SR$exonBound5p=exonBound5p
      myFusionFinal.SR$exonID3p=exonID3p
      myFusionFinal.SR$exonBound3p=exonBound3p
    }
  }

##### do similarly for Mapped reads

 #make columns first
  myFusionFinal.MR$exonID5p= myFusionFinal.MR$brpos5.start
  myFusionFinal.MR$exonBound5p= myFusionFinal.MR$brpos5.start
  myFusionFinal.MR$exonID3p= myFusionFinal.MR$brpos5.start
  myFusionFinal.MR$exonBound3p= myFusionFinal.MR$brpos5.start

  if (nrow(myFusionFinal.MR) > 0){
    myFusionFinal.MR$exonID5p=-1
    myFusionFinal.MR$exonBound5p=-1
    myFusionFinal.MR$exonID3p=-1
    myFusionFinal.MR$exonBound3p=-1
    if (nrow(myFusionFinal.MR)>0){
      exonInfo=select(anntxdb, keys=unique(c(as.character(myFusionFinal.MR$gene5p),as.character(myFusionFinal.MR$gene3p))), columns=c("GENEID","EXONID","EXONSTART","EXONEND"), keytype = "GENEID")
      exonBound5p=exonBound3p=NULL
      exonID5p=exonID3p=NULL
      for (i in 1:nrow(myFusionFinal.MR)){
        x=myFusionFinal.MR[i,]  
        br5=x$brpos5.start
        g5=as.character(x$gene5p)
        e5=exonInfo[exonInfo$GENEID == g5,]
        if (x$strand5p=="+"){
          absdif=abs(br5-e5$EXONEND)
          k=which.min(absdif)
          exonID5p=c(exonID5p,e5$EXONID[k])
          exonBound5p=c(exonBound5p,e5$EXONEND[k])
        }
        if (x$strand5p=="-"){
          absdif=abs(br5-e5$EXONSTART)
          k=which.min(absdif)
          exonID5p=c(exonID5p,e5$EXONID[k])
          exonBound5p=c(exonBound5p,e5$EXONSTART[k])    
        }
        br3=x$brpos3.start
        g3=as.character(x$gene3p)
        e3=exonInfo[exonInfo$GENEID == g3,]
        if (x$strand3p=="+"){
          absdif=abs(br3-e3$EXONSTART)
          k=which.min(absdif)
          exonID3p=c(exonID3p,e3$EXONID[k])
          exonBound3p=c(exonBound3p,e3$EXONSTART[k])
        }
        if (x$strand3p=="-"){
          absdif=abs(br3-e3$EXONEND)
          k=which.min(absdif)
          exonID3p=c(exonID3p,e3$EXONID[k])
          exonBound3p=c(exonBound3p,e3$EXONEND[k])
        }
      }
      myFusionFinal.MR$exonID5p=exonID5p
      myFusionFinal.MR$exonBound5p=exonBound5p
      myFusionFinal.MR$exonID3p=exonID3p
      myFusionFinal.MR$exonBound3p=exonBound3p
    }
  }

  ##### merging two fusion lists to be one
  selectedColumns=c("chrom5p","brpos5.start","brpos5.end","chrom3p","exonBound5p","exonID5p","brpos3.start","brpos3.end","exonBound3p","exonID3p","fusionName","score","strand5p","strand3p","supportRead")
  myFusionFinal=myFusionFinal.SR[,selectedColumns]
  if (nrow(myFusionFinal.MR)>0){
    newID=myFusionFinal.MR$fusionName %in% myFusionFinal$fusionName
    myFusionFinal=rbind(myFusionFinal,myFusionFinal.MR[!newID,selectedColumns])
  }
  
  #myFusionFinal=myFusionFinal[!duplicated(myFusionFinal$fusionName),]
  dim(myFusionFinal)
  myFusionFinal$name12=myFusionFinal$fusionName
  myFusionFinal$name21=sapply(as.character(myFusionFinal$fusionName),function(x) paste(rev(unlist(strsplit(x,"-"))),collapse ="-"))
  
  myFusionFinal$gene5=sapply(as.character(myFusionFinal$fusionName),function(x) unlist(strsplit(x,"-"))[1])
  myFusionFinal$gene3=sapply(as.character(myFusionFinal$fusionName),function(x) unlist(strsplit(x,"-"))[2])
  
  #get extra information of genes
  myFusionFinal$ribSub5=myFusionFinal$gene5;myFusionFinal$ribSub3=myFusionFinal$gene5;myFusionFinal$mitoTrans5=myFusionFinal$gene5;myFusionFinal$mitoTrans3=myFusionFinal$gene5;myFusionFinal$ribonupro5=myFusionFinal$gene5;myFusionFinal$ribonupro3=myFusionFinal$gene5;

  if (nrow(myFusionFinal)>0){
    myFusionFinal$ribSub5=NA;myFusionFinal$ribSub3=NA;myFusionFinal$mitoTrans5=NA;myFusionFinal$mitoTrans3=NA;myFusionFinal$ribonupro5=NA;myFusionFinal$ribonupro3=NA;
    #if they are from cytosolic ribosomal subunit
    if (exists("ribSubunitDb")){
      myFusionFinal$ribSub5=match(as.character(myFusionFinal$gene5),as.character(ribSubunitDb$ensembl_gene_id))
      myFusionFinal$ribSub3=match(as.character(myFusionFinal$gene3),as.character(ribSubunitDb$ensembl_gene_id))
    }
    #if they are from mitochondrial transclation
    if (exists("mitoTransDb")){
      myFusionFinal$mitoTrans5=match(as.character(myFusionFinal$gene5),as.character(mitoTransDb$ensembl_gene_id))
      myFusionFinal$mitoTrans3=match(as.character(myFusionFinal$gene3),as.character(mitoTransDb$ensembl_gene_id))
    }  
    #if they are from ribonucleoprotein complex
    if (exists("ribonuproDb")){
      myFusionFinal$ribonupro5=match(as.character(myFusionFinal$gene5),as.character(ribonuproDb$ensembl_gene_id))
      myFusionFinal$ribonupro3=match(as.character(myFusionFinal$gene3),as.character(ribonuproDb$ensembl_gene_id))
    }
  }
  
  #order by score
  myFusionFinal=myFusionFinal[order(myFusionFinal$score,decreasing = TRUE),]
  
  if (nrow(myFusionFinal)==0){
    return(list(myFusionFinal=myFusionFinal,myFusionFinal.SR=myFusionFinal.SR,myFusionFinal.MR=myFusionFinal.MR))
  }
  
  myFusionFinal$MR=myFusionFinal$SR=0
  #myFusionFinal$SR[which(!is.na(match(as.character(myFusionFinal$fusionName),as.character(myFusionFinal.SR$fusionName))))]=1
  #myFusionFinal$MR[which(!is.na(match(as.character(myFusionFinal$fusionName),as.character(myFusionFinal.MR$fusionName))))]=1  
  matchID=match(as.character(myFusionFinal$fusionName),as.character(myFusionFinal.SR$fusionName))
  pick=which(!is.na(matchID))
  myFusionFinal$SR[pick]=myFusionFinal.SR$supportCount[matchID[pick]]
  matchID=match(as.character(myFusionFinal$fusionName),as.character(myFusionFinal.MR$fusionName))
  pick=which(!is.na(matchID))
  myFusionFinal$MR[pick]=myFusionFinal.MR$supportCount[matchID[pick]]

  
  return(list(myFusionFinal=myFusionFinal,myFusionFinal.SR=myFusionFinal.SR,myFusionFinal.MR=myFusionFinal.MR))
}




