#This script produces a table of unique alleles with their abundance accross time points,including information on the locus where it falls, the coordinates of the locus and its orientation. It requires the pipeline main_miseq.sh to be executed first.

#FIRST, LOAD THE INSERTION ALLELE TABLES
S0<-read.table('R0_results.txt', sep="\t", header=FALSE, as.is=TRUE)
S1<-read.table('R1_results.txt', sep="\t", header=FALSE, as.is=TRUE)
S3<-read.table('R3_results.txt', sep="\t", header=FALSE, as.is=TRUE)
S5<-read.table('R5_results.txt', sep="\t", header=FALSE, as.is=TRUE)
S8<-read.table('R8_results.txt', sep="\t", header=FALSE, as.is=TRUE)

#AGGREGATE DATA FROM DIFFERENT TIME POINTS
tot<-rbind(S0[,1:3], S1[,1:3], S3[,1:3], S5[,1:3], S8[,1:3])
unlist<-unique(tot[,1:2])		#compile a list of unique insertion alleles
table<-cbind(unlist, 0,0,0,0,0,0)	#create empty table to record abundance

#SUM ABUNDANCES OF UNIQUE ALLELES WITHIN EACH TIMEPOINT
for (i in 1:dim(unlist)[1]) {
    if (i%%10000==0) {print(c(i,dim(unlist)))}
    if (unlist[i,2]%in%S0$V2) {
        table[i,3] <- sum(S0[S0$V2==unlist[i,2],3])
    } else {table[i,3] <- 0}

    if (unlist[i,2]%in%S1$V2) {
        table[i,4] <- sum(S1[S1$V2==unlist[i,2],3])
    } else {table[i,4] <- 0}

    if (unlist[i,2]%in%S3$V2) {
        table[i,5] <- sum(S3[S3$V2==unlist[i,2],3])
    } else {table[i,5] <- 0}

    if (unlist[i,2]%in%S5$V2) {
        table[i,6] <- sum(S5[S5$V2==unlist[i,2],3])
    } else {table[i,6] <- 0}

    if (unlist[i,2]%in%S8$V2) {
        table[i,7] <- sum(S8[S8$V2==unlist[i,2],3])
    } else {table[i,7] <- 0}

}

#ADD INFORMATION ON LOCUS
index<-read.table('parsed_R606genoscope_I.txt', sep="\t", header=FALSE, as.is=TRUE)
table<-cbind(table, 0,0,0)
colnames(table)<-c('name', 'alle', 't0', 't1', 't3', 't5', 't8', 'strand', 'srt', 'end')
rownames(table)<-seq(1,dim(unlist)[1],1)

for (i in 1:dim(unlist)[1]) {

    table[i,8:10]<-index[index[,1]==table[i,1],c(4,2,3)]	#record locus information 
    if (i%%1000==0) {print(i)}
}

table<-table[order(table$alle),]					#sort the table for downstream analyses
write.table(table, file = "Rtable(fuzz).txt", sep = "\t")		#write file
source('nplicates_updater.R') 					#to handle the issues arising from n-plicated ORFs
