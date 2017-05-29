############################################################
#####extract biological features of fusion gene candidates in mapped read pipeline and do some filters
#####input: geneParalog, hgncName, ribSubunitDb, mitoTransDb, ribonuproDb ... from global enviroment loaded from txAnnofile beforeward
doBiologicalFilter<-function(fgeList, chromRef=paste("",c(1:22,"X","Y"),sep=""),onlyProteinCodingGenes=TRUE, doFilter=TRUE){
  myFusionFinal=fgeList
  
  if (onlyProteinCodingGenes){
    cat("\n Keep only protein-coding genes")
    keepID=which(myFusionFinal$geneType1=="protein_coding" & myFusionFinal$geneType2=="protein_coding")
    myFusionFinal=myFusionFinal[keepID,]
  }
    
  rmID=unlist(lapply(c(1:nrow(myFusionFinal)), function(x){
    par1=c(as.character(myFusionFinal$gene1[x]),geneParalog[which(geneParalog[,1]==as.character(myFusionFinal$gene1[x])),2])
    par2=c(as.character(myFusionFinal$gene2[x]),geneParalog[which(geneParalog[,1]==as.character(myFusionFinal$gene2[x])),2])
    if (length(par1)>0 & length(par2)>0) return(length(intersect(par1,par2))>0)
    return(FALSE)
  }))
  myFusionFinal$paralog=rep(0,nrow(myFusionFinal))
  myFusionFinal$paralog[rmID]=1

  if (doFilter){
    cat("\n Eliminate the fusion between genes and their paralogs")
    myFusionFinal=myFusionFinal[myFusionFinal$paralog==0,]
  }
  
  gene5=as.character(myFusionFinal$gene1)
  res <- hgncName[match(gene5,as.character(hgncName$ensembl_gene_id)) ,]
  res=res[res$chromosome_name%in%chromRef,]
  matchID=match(gene5,res$ensembl_gene_id)
  myFusionFinal$gene1_ucsc=res$hgnc_symbol[matchID]
  gene3=as.character(myFusionFinal$gene2)
  res <- hgncName[match(gene3,as.character(hgncName$ensembl_gene_id)) ,]
  res=res[res$chromosome_name%in%chromRef,]
  matchID=match(gene3,res$ensembl_gene_id)
  myFusionFinal$gene2_ucsc=res$hgnc_symbol[matchID]
   
  #filter a fusion made by a gene and its read-through genes
  keepID=unlist(lapply(c(1:nrow(myFusionFinal)), function(x){
    res=c(unlist(strsplit(myFusionFinal$gene1_ucsc[x],"-")),unlist(strsplit(myFusionFinal$gene2_ucsc[x],"-")))
    return((length(res) <= 2))
  }))
  
  readThrough=rep(1,nrow(myFusionFinal))
  if (length(keepID)>0) readThrough[keepID]=0
  myFusionFinal$readThrough=readThrough
  if (doFilter){
    cat("\n Eliminate the fusion between single gene and read-through gene")
    myFusionFinal=myFusionFinal[myFusionFinal$readThrough==0,]
    
  }

  #compute again multiple duplicated genes
  res=computeDupGene(myFusionFinal,dupGene.thres=-1)
  colnames(res)=c("dupGene1_f1","dupGene2_f1")
  myFusionFinal=cbind(myFusionFinal,res)
  
  return(myFusionFinal)
}