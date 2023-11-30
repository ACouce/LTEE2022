#Here is the algorithm that pools alleles at the level of sub-genic regions

# .. make a column to label the different alleles (i.e, alias)
alls$fitness<-alls$name
colnames(alls)[colnames(alls) == 'fitness'] <- 'alias'

# .. reorganize columns for convenience downstream
alls<-alls[,c(1,11,2,3,4,5,6,7,8,9,10)]

# .. the algorithm
if (dim(alls)[1]>=1) {	# do we have at least 1 allele?

	alls$sum<-rowSums(alls[,4:8]) 			# make a column with sum of observations across time
	alls<-subset(alls[alls$sum>=3,])		 	# pick up all alleles abundant enough (threshold from tests with neutrals)

	if (dim(alls)[1]>=1) {	# do we have at least 1 allele?
		pool<-alls
		pool$pos<-(pool$alle-pool$srt)/(pool$end-pool$srt)	# record position within unit
		pool<-pool[,-(10:11)]					# remove columns with start and end of sub-genic region information
		pool$alias<-orf						# record name of sub-genic region
		row.names(pool)<-NULL					# polish data frame
		pool$abc<-dim(alls)[1]					# record number of different alleles
		pool<-pool[,c(1,2,3,12,4,5,6,7,8,10,11,9)]		# reorganize columns for convenience downstream
		pool1<-pool[1,]						# prepare to add up all alleles
		pool1[1,5:10]<-colSums(pool[,5:10])			# the addition	
		list<-rbind(list, pool1)					# update recyclabe list
		listtot<-rbind(listtot, list)				# update global list	
	}
	
} else {	# to handle sub-genic regions without observations

	list[1:2]<-orf

}
