# .. this script generates simulated DFEs by randomly sampling a quantity of loci equal to the observed number of beneficial mutations in the real DFE, providing a control for neutral expectation.

# .. set number of simulations
replicates<-100

# .. prepare the pool of neutral and deleterious mutations to sample with (for the ancestor)
Rtable<-read.table(file = "Rfitted_log.txt", sep="\t", header=TRUE, as.is=TRUE, comment.char = "")
Rtable<-Rtable[!is.na(Rtable$sterr1),]
quasiN<-Rtable[which(Rtable$fitted1>-0.3 & Rtable$fitted1<=0), c('ORF','fitted1')] 
colnames(quasiN)<-NULL

# .. record the number of beneficial mutations observed in the ancestor (to ensure proper control)
size<-nrow(read.table(file = 'Rtoanal.dat', sep = "\t", header = TRUE))

# .. create empty data frame to store simulated data
Ntabs<-as.data.frame(matrix(0, size*replicates, 3))
colnames(Ntabs) <- c('name', 'fitted1', 'rep')

# .. the sampling algorithm
c<-0
for (r in 1:replicates) {	
	sample<-sample(1:nrow(quasiN), size, replace = FALSE, prob = NULL)	# sampling without replacement
	str<-(r-1)*size+1							# trick to write the blocks of data in consecutive rows
	end<-r*size
	Ntabs[str:end,c('name','fitted1')]<-quasiN[sample,]			# store sampled fitness effects
	Ntabs[str:end,c('rep')]<-r						# add label for replica number
}

# .. write simulated data to be used later by 'representativeness_sim.R'
write.table(Ntabs, 'Rcontrol.dat', row.names = FALSE, sep = "\t",  col.names = TRUE)


# .. prepare the pool of neutral and deleterious mutations to sample with (for 2K)
Rtable<-read.table(file = "2Kfitted_log.txt", sep="\t", header=TRUE, as.is=TRUE, comment.char = "")
Rtable<-Rtable[!is.na(Rtable$sterr1),]
quasiN<-Rtable[which(Rtable$fitted1>-0.3 & Rtable$fitted1<=0), c('ORF','fitted1')] 
colnames(quasiN)<-NULL

# .. record the number of beneficial mutations observed in 2K (to ensure proper control)
size<-nrow(read.table(file = 'Ttoanal.dat', sep = "\t", header = TRUE))

# .. create empty data frame to store simulated data
Ntabs<-as.data.frame(matrix(0, size*replicates, 3))
colnames(Ntabs) <- c('name', 'fitted1', 'rep')

# .. the sampling algorithm
c<-0
for (r in 1:replicates) {	
	sample<-sample(1:nrow(quasiN), size, replace = FALSE, prob = NULL)	# sampling without replacement
	str<-(r-1)*size+1							# trick to write the blocks of data in consecutive rows
	end<-r*size
	Ntabs[str:end,c('name','fitted1')]<-quasiN[sample,]			# store sampled fitness effects
	Ntabs[str:end,c('rep')]<-r						# add label for replica number
}
	
# .. write simulated data to be used later by 'representativeness_sim.R'	
write.table(Ntabs, 'Tcontrol.dat', row.names = FALSE, sep = "\t",  col.names = TRUE)


# .. prepare the pool of neutral and deleterious mutations to sample with (for 15K)
Rtable<-read.table(file = "15Kfitted_log.txt", sep="\t", header=TRUE, as.is=TRUE, comment.char = "")
Rtable<-Rtable[!is.na(Rtable$sterr1),]
quasiN<-Rtable[which(Rtable$fitted1>-0.3 & Rtable$fitted1<=0), c('ORF','fitted1')] 
colnames(quasiN)<-NULL

# .. record the number of beneficial mutations observed in 15K (to ensure proper control)
size<-nrow(read.table(file = 'Ftoanal.dat', sep = "\t", header = TRUE))

#create empty data frame
Ntabs<-as.data.frame(matrix(0, size*replicates, 3))
colnames(Ntabs) <- c('name', 'fitted1', 'rep')

# .. the sampling algorithm
c<-0
for (r in 1:replicates) {	
	sample<-sample(1:nrow(quasiN), size, replace = FALSE, prob = NULL)	# sampling without replacement
	str<-(r-1)*size+1							# trick to write the blocks of data in consecutive rows
	end<-r*size
	Ntabs[str:end,c('name','fitted1')]<-quasiN[sample,]			# store sampled fitness effects
	Ntabs[str:end,c('rep')]<-r						# add label for replica number
}

# .. write simulated data to be used later by 'representativeness_sim.R'	
write.table(Ntabs, 'Fcontrol.dat', row.names = FALSE, sep = "\t",  col.names = TRUE)
