#This script produces some statistics and prepares working objects


# .. read the data after pooling
ttable<-read.table(file = filename2, sep="\t", header=TRUE, as.is=TRUE, comment.char = "") #The default is comment.char = "#"!!

# .. prepare working objects
alleles<-length(unique(ttable$alle))
rawcounts<-matrix(0,alleles,5)	
rawcounts<-ttable[,5:9]
rawstat<-matrix(0,5,6)

# .. loop to estimate basic statistics of the raw data
for (i in 1:5) {
	rawstat[i,1]<-sum(rawcounts[,i], na.rm='T') 
	rawstat[i,2]<-sum(length(unique(ttable[ttable[i+2]>0,]$ORF)))
 	rawstat[i,3]<-sum(length(unique(ttable[ttable[i+2]>0,]$alle)))
 	rawstat[i,4]<-sum(rawcounts[,i]==0, na.rm='T')
	rawstat[i,5]<-median(rawcounts[,i], na.rm='T')
    	rawstat[i,6]<-mean(rawcounts[,i], na.rm='T')
}

colnames(rawstat)<-c("total","loci","alleles","zeros","median","average")
print(rawstat)
