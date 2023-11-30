#this master script compares different datasets and normalizes their coverages via resampling according to median number of insertions in the ancestor. The actual sample files used in the analyses are provided.

# .. remove any previous content in the session
rm(list=ls())

# .. read the data 
Rtable<-read.table(file = 'Rtable_updated(fuzz).txt', sep="\t", header=TRUE, as.is=TRUE, comment.char = "") #The default is comment.char = "#"!!)
Ttable<-read.table(file = 'Ttable_updated(fuzz).txt', sep="\t", header=TRUE, as.is=TRUE, comment.char = "") #The default is comment.char = "#"!!)
Ftable<-read.table(file = 'Ftable_updated(fuzz).txt', sep="\t", header=TRUE, as.is=TRUE, comment.char = "") #The default is comment.char = "#"!!)

# .. create matrix to loop and estimate basic statistics of the data
rawstat<-matrix(0,5,3)
for (i in 1:5) {
	rawstat[i,1]<-sum(Rtable[,i+2], na.rm='T')
	rawstat[i,2]<-sum(Ttable[,i+2], na.rm='T') 
	rawstat[i,3]<-sum(Ftable[,i+2], na.rm='T') 
}

rawstat[,2]<-floor(rawstat[,2]*c(1.33,0.73,0.8,0.44,0.56)) # ratio of actual median number of inserstions per locus 2K vs Anc
rawstat[,3]<-floor(rawstat[,3]*c(1,0.69,0.62,0.44,0.6))	# ratio of actual median number of inserstions per locus 15K vs Anc


# .. resampling to normalize algorithm (10 replicates)
for (t in 3:7) {
	for (r in 1:10) {	
		sampling<-sample(Ttable_ori$alle, rawstat[t-2,2], replace = TRUE, prob=Ttable_ori[,t])
		tab<-table(sampling)
		alleles<-as.numeric(names(tab))
		Ttable[Ttable$alle%in%alleles,t]<-Ttable[Ttable$alle%in%alleles,t]+as.numeric(tab)
		}
	}
	
	# .. resampling to normalize algorithm (10 replicates)
for (t in 3:7) {
	for (r in 1:10) {	
		sampling<-sample(Ftable_ori$alle,  rawstat[t-2,3], replace = TRUE, prob=Ftable_ori[,t])
		tab<-table(sampling)
		alleles<-as.numeric(names(tab))
		Ftable[Ftable$alle%in%alleles,t]<-Ftable[Ftable$alle%in%alleles,t]+as.numeric(tab)
		}
	}

# .. divide by the number of replicates
Ttable[,3:7]<-floor(Ttable[,3:7]/10)
Ftable[,3:7]<-floor(Ftable[,3:7]/10)
	
# .. save results
write.table(Rtable, file = "Rtable_updated(norm).txt", sep = "\t")	 # add this one for convenience downstream
write.table(Ttable, file = "Ttable_updated(norm).txt", sep = "\t")
write.table(Ftable, file = "Ftable_updated(norm).txt", sep = "\t")
