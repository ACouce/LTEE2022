#this script merges the lists of UMI-corrected alleles, and prepares a table with useful information for each one (e.g. names, time series, strand, genomic coordinates)

# .. load genomic coordinates
index<-read.table('parsed_R606genoscope_I.txt', sep="\t", header=FALSE, as.is=TRUE)

# .. load lists of UMI-corrected alleles
J1<-read.table('ara+2_J1_allcount.txt', header=FALSE, as.is=TRUE)
J2<-read.table('ara+2_J2_allcount.txt', header=FALSE, as.is=TRUE)
J3<-read.table('ara+2_J3_allcount.txt', header=FALSE, as.is=TRUE)
J4<-read.table('ara+2_J4_allcount.txt', header=FALSE, as.is=TRUE)
J5<-read.table('ara+2_J5_allcount.txt', header=FALSE, as.is=TRUE)

# .. aggregate data from different time points
tot<-rbind(J1[,1:3], J2[,1:3], J3[,1:3], J4[,1:3], J5[,1:3])
unlist<-unique(tot[,1:2])				#compile a list of unique insertion alleles
unlist<-unlist[order(unlist$V2),]
table<-cbind(unlist, 0,0,0,0,0,0,0,0)			#create empty table to record values

uniques<-dim(unlist)[1]

for (i in 1:uniques) {
    if (i%%10000==0) {print(c(i,dim(unlist)))}

    if (unlist[i,2]%in%J1$V2) {table[i,3] <- J1[J1$V2==unlist[i,2],4]}	# record abundance values
    if (unlist[i,2]%in%J2$V2) {table[i,4] <- J2[J2$V2==unlist[i,2],4]}
    if (unlist[i,2]%in%J3$V2) {table[i,5] <- J3[J3$V2==unlist[i,2],4]}
    if (unlist[i,2]%in%J4$V2) {table[i,6] <- J4[J4$V2==unlist[i,2],4]}
    if (unlist[i,2]%in%J5$V2) {table[i,7] <- J5[J5$V2==unlist[i,2],4]}
    
    table[i,8:10]<-index[index[,1]==table[i,1],c(4,2,3)]	#add information on strand, start and end
}        

colnames(table)<-c('name', 'alle', 't1', 't2', 't3', 't4', 't5', 'strand', 'srt', 'end')
rownames(table)<-seq(1,uniques,1)

table<-table[order(table$alle),]				# sort table for downstream convenience

# .. write resulting table
write.table(table, file = "2Ktable(fuzz).txt", sep = "\t")
