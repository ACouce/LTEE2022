# .. This script filters out low-quality fits and presumed hitchhikers, and then constructs the count tables required for plotting the DFE.

# .. define a useful function
'%!in%' <- function(x,y)!('%in%'(x,y))

# .. load the specific dataset to work with
ttable<-read.table(file = filename2, sep="\t", header=TRUE, as.is=TRUE, comment.char = "")

# .. remove NAs, loci with very low initial counts, less than two independent alleles, or large regression standard errors

ttable[is.na(ttable$pval),c(4:24)]<--100
ttable<-ttable[!is.na(ttable$fitted1),]
ttable[(ttable$t2<=3 | ttable$t5<=15) & ttable$fitted1>0,]$sterr1 <-NA
ttable<-ttable[!is.na(ttable$sterr1),]
ttable[ttable$sterr1>0.012 & ttable$fitted1>0 ,]$sterr1 <-NA
ttable[ttable$abn<2,]$sterr1 <-NA
ttable<-ttable[!is.na(ttable$sterr1),]

# .. store clean data table for further reference
write.table(ttable, file = filename_fil, col.names = TRUE, sep="\t")

# .. create objects for constructing a count table of fitness classes (required for plotting the DFE)
r_fitted_count<-matrix(0,nrow = 60, ncol = 4)
orgn=-0.33
step=0.015
for (i in 2:60) {	# loop accross fitness classes defined by orgn and step
	low=round(orgn+step*i, digits = 3)
	hgh=round(low+step, digits = 3)
	r_fitted_count[i,1]<-round(mean(c(low,hgh)), digits = 3)
	r_fitted_count[i,2]<-sum(ttable$fitted>=low, na.rm=TRUE)-sum(ttable$fitted>hgh, na.rm=TRUE)
	r_fitted_count[i,3]<-sum(ttable$fitted1>=low, na.rm=TRUE)-sum(ttable$fitted1>hgh, na.rm=TRUE)
	r_fitted_count[i,4]<-r_fitted_count[i-1,4]+r_fitted_count[i,3]
		}
		
# .. store count table of fitness classes
write.table(r_fitted_count, paste(label,'hist.dat', sep="", collapse=NULL), col.names = FALSE, row.names = FALSE)

# .. print best 10 alleles
best<-ttable[ttable$fitted1>0.02,]
chosen<-best[!is.na(best$t1),c(2,4:11, 14, 17,20,23)]
chosen<-chosen[order(chosen$fitted1, decreasing = TRUE),]
print(head(chosen, 10))
