############################################################
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
  #keep only one presentative for each fusion gene
  myFusionFinal.SR=myFusionFinal.SR[!duplicated(myFusionFinal.SR$name12),]
  
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
  
  
  
  #create common features
  myFusionFinal.SR$chrom5p=myFusionFinal.SR$chrom1
  myFusionFinal.SR$brpos5.start=myFusionFinal.SR$brchposEx5
  myFusionFinal.SR$brpos5.end=myFusionFinal.SR$brchposEx5+1
  myFusionFinal.SR$strand5p=myFusionFinal.SR$strand1
  myFusionFinal.SR$chrom3p=myFusionFinal.SR$chrom2
  myFusionFinal.SR$brpos3.start=myFusionFinal.SR$brchposEx3
  myFusionFinal.SR$brpos3.end=myFusionFinal.SR$brchposEx3+1
  myFusionFinal.SR$strand3p=myFusionFinal.SR$strand2

  ##### merging two fusion lists to be one
  myFusionFinal=myFusionFinal.SR[,c("chrom5p","brpos5.start","brpos5.end","chrom3p","brpos3.start","brpos3.end","fusionName","score","strand5p","strand3p","supportRead")]
  if (nrow(myFusionFinal.MR)>0) myFusionFinal=rbind(myFusionFinal,myFusionFinal.MR[,c("chrom5p","brpos5.start","brpos5.end","chrom3p","brpos3.start","brpos3.end","fusionName","score","strand5p","strand3p","supportRead")])
  
  myFusionFinal=myFusionFinal[!duplicated(myFusionFinal$fusionName),]
  dim(myFusionFinal)
  myFusionFinal$name12=myFusionFinal$fusionName
  myFusionFinal$name21=sapply(as.character(myFusionFinal$fusionName),function(x) paste(rev(unlist(strsplit(x,"-"))),collapse ="-"))
  
  myFusionFinal$gene5=sapply(as.character(myFusionFinal$fusionName),function(x) unlist(strsplit(x,"-"))[1])
  myFusionFinal$gene3=sapply(as.character(myFusionFinal$fusionName),function(x) unlist(strsplit(x,"-"))[2])
  
  
  #if they are from cytosolic ribosomal subunit
  myFusionFinal$ribSub5=match(as.character(myFusionFinal$gene5),as.character(ribSubunitDb$ensembl_gene_id))
  myFusionFinal$ribSub3=match(as.character(myFusionFinal$gene3),as.character(ribSubunitDb$ensembl_gene_id))
  #if they are from mitochondrial transclation
  myFusionFinal$mitoTrans5=match(as.character(myFusionFinal$gene5),as.character(mitoTransDb$ensembl_gene_id))
  myFusionFinal$mitoTrans3=match(as.character(myFusionFinal$gene3),as.character(mitoTransDb$ensembl_gene_id))
  #if they are from ribonucleoprotein complex
  myFusionFinal$ribonupro5=match(as.character(myFusionFinal$gene5),as.character(ribonuproDb$ensembl_gene_id))
  myFusionFinal$ribonupro3=match(as.character(myFusionFinal$gene3),as.character(ribonuproDb$ensembl_gene_id))
  
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



