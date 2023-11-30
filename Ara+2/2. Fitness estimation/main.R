#This master script handles the pipeline for estimating fitness from insertion count time series.

# .. remove any previous content in the session
rm(list=ls())

# .. start the stopwatch
clock <- Sys.time()

# .. matrix to help calling and naming the different data files along the loop
labels<-matrix(data=c('R','Ap2_2K','Ap2_15K'), nrow=1) 

for (sour in 1:3) {	# loop through the time series data files
	
	label<-labels[sour]
	filename1<-paste(label,"table_updated(fuzz).txt",sep="")	# to read insertion count time series
	filename2<-paste(label,"list_total_filtered.txt",sep="")	# to store pooled insertion counts
	filename3<-paste(label,"fitted_log.txt",sep="")		# to stores fitness estimates from pooled insertion frequencies 

	source('allele_pooling.R')	# produces time series files with the pooled data at the level of sub-genic regions (filename2)
	source('stats.R')		# produces some statistics and prepares working objects
	source('fitting.R')		# produces files with the estimated fitnesses (filename3)
	
	}

# .. remove intermediate files
file.remove(list.files(pattern = "*filtered.txt"))

# .. stop the stopwatch
print(Sys.time()-clock)
