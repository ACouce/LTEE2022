#This script produces a table of unique alleles with their abundance accross time points,including information on the locus where it falls, the coordinates of the locus and its orientation.
#This script requires the pipeline main.sh to be executed first. In the confuguiration here presented, the script will merge reads from independent runs of the same samples (e.g., Jan 2017 and Sep 2017).

#FIRST, LOAD THE INSERTION ALLELE TABLES
Jan1<-read.table('Jan17_Ap2_2K1_results.txt', sep="\t", header=FALSE, as.is=TRUE)
Jan2<-read.table('Jan17_Ap2_2K2_results.txt', sep="\t", header=FALSE, as.is=TRUE)
Jan3<-read.table('Jan17_Ap2_2K3_results.txt', sep="\t", header=FALSE, as.is=TRUE)
Jan4<-read.table('Jan17_Ap2_2K4_results.txt', sep="\t", header=FALSE, as.is=TRUE)
Jan5<-read.table('Jan17_Ap2_2K5_results.txt', sep="\t", header=FALSE, as.is=TRUE)

Sep1<-read.table('Sep17_Ap2_2K1_results.txt', sep="\t", header=FALSE, as.is=TRUE)
Sep2<-read.table('Sep17_Ap2_2K2_results.txt', sep="\t", header=FALSE, as.is=TRUE)
Sep3<-read.table('Sep17_Ap2_2K3_results.txt', sep="\t", header=FALSE, as.is=TRUE)
Sep4<-read.table('Sep17_Ap2_2K4_results.txt', sep="\t", header=FALSE, as.is=TRUE)
Sep5<-read.table('Sep17_Ap2_2K5_results.txt', sep="\t", header=FALSE, as.is=TRUE)

#MERGE DATA FROM SAME SAMPLES
merged1<-rbind(Jan1[,1:3], Sep1[,1:3])
merged2<-rbind(Jan2[,1:3], Sep2[,1:3))
merged3<-rbind(Jan3[,1:3], Sep3[,1:3])
merged4<-rbind(Jan4[,1:3], Sep4[,1:3])
merged5<-rbind(Jan5[,1:3], Sep5[,1:3])

#AGGREGATE DATA FROM DIFFERENT TIME POINTS
tot<-rbind(merged1[,1:3], merged2[,1:3], merged3[,1:3], merged4[,1:3], merged5[,1:3])

#PREPARE A TABLE TO RECORD ABUNDANCE ACCROSS TIME POINTS
unlist<-unique(tot[,1:2])		#compile a list of unique insertion alleles
table<-cbind(unlist, 0,0,0,0,0)		#create empty table to record abundance

#SUM ABUNDANCES OF UNIQUE ALLELES WITHIN EACH TIMEPOINT
for (i in 1:dim(unlist)[1]) {
    if (i%%10000==0) {print(c(i,dim(unlist)))}

    if (unlist[i,2]%in%merged1$V2) {
        table[i,4] <- sum(merged1[merged1$V2==unlist[i,2],3])
    } else {table[i,4] <- 0}

    if (unlist[i,2]%in%merged2$V2) {
        table[i,5] <- sum(merged2[merged2$V2==unlist[i,2],3])
    } else {table[i,5] <- 0}

    if (unlist[i,2]%in%merged3$V2) {
        table[i,6] <- sum(merged3[merged3$V2==unlist[i,2],3])
    } else {table[i,6] <- 0}

    if (unlist[i,2]%in%merged4$V2) {
        table[i,7] <- sum(merged4[merged4$V2==unlist[i,2],3])
    } else {table[i,7] <- 0}

    if (unlist[i,2]%in%merged5$V2) {
        table[i,8] <- sum(merged5[merged5$V2==unlist[i,2],3])
    } else {table[i,8] <- 0}
}

#ADD INFORMATION ON LOCUS
index<-read.table('parsed_R606genoscope_IU.txt', sep="\t", header=FALSE, as.is=TRUE)
table<-cbind(table, 0,0,0)		#create empty columns to store locus information 
colnames(table)<-c('name', 'alle', 't1', 't2', 't3', 't4', 't5', 'strand', 'srt', 'end')
rownames(table)<-seq(1,dim(unlist)[1],1)

for (i in 1:dim(unlist)[1]) {
    table[i,8:10]<-index[index[,1]==table[i,1],c(4,2,3)]	#record locus information 
    if (i%%1000==0) {print(i)}
}

table<-table[order(table$alle),]		#sort the table for downstream analyses
write.table(table, file = "Ap2_2Ktable_updated(fuzz).txt", sep = "\t")	#write file
