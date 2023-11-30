#this script takes the 'merged.txt' files (with format: locus position UMI) and produces unique lists of alleles (whith format: locus raw reads UMI-corrected reads).

# .. start the stopwatch
clock <- Sys.time()

suffix<-'ara2-2K_'

for (day in c('J1','J2','J3','J4','J5')) {

	x<-read.table(paste0(suffix,day,'_merged.txt'))
	tab<-table(x$V3)

	df<-as.data.frame(matrix(NA,dim(tab),4))
	for (n in 1:dim(tab)) {

		df[n,1]<-as.character(x[x$V3%in%as.numeric(names(tab[n])),][1,2])
		df[n,2]<-as.numeric(names(tab[n]))
		df[n,3]<-as.numeric(tab[n])
		
		# only consider umis if appropriate:
		if (df[n,3]<=1) {
			df[n,4]<-df[n,3]
		} else {
			df[n,4]<-length(unique(x[x$V3==df[n,2],]$V4))
		}
		
		if (n%%5000==0) {print(n)}	
	}
	write.table(df,paste0(suffix,day,'_allcount.txt'), quote = FALSE, col.names = FALSE, row.names= FALSE)
}

# .. stop the stopwatch
print(Sys.time()-clock)
