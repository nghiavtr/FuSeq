##### Date: 08 Nov 2018
##### Contact: Trung Nghia Vu (nghiavtr@gmail.com)
##### This script is to remove version number in the transcript name from ensembl annotation

args = commandArgs(trailingOnly=TRUE)
cdnaFn=args[1]
outFn=args[2]

#read fasta file
con <- file(cdnaFn, "r", blocking = FALSE)
mydata=readLines(con)
close(con)

txheaderID=grep(">",mydata)
#remove version in 
txName=sapply(txheaderID, function(i){
	x=mydata[i]
	y=unlist(strsplit(x," "))
  	t=substring(y[1],2)
  	t2=unlist(strsplit(t,"\\."))[1]
  	x1=gsub(t,t2,x)
  	return(x1)
})

#replace with new txname
mydata[txheaderID]=txName
#export to file
write(mydata,outFn)


