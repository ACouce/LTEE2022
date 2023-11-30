#This script loops through the three datasets cleansing out low-quality fits and presumed hitchhikers, subsequently constructing the count tables required for plotting the DFE.

# .. load a library useful for some data manipulations (e.g., between()) 
library("dplyr")

# .. matrices to help calling and naming the different data files along the loop
files<-matrix(data=c("Rfitted_log.txt", "2Kfitted_log.txt", "15Kfitted_log.txt"),nrow=1)		#the input files
files_fil<-matrix(data=c("Rfitted_fil.txt", "2Kfitted_fil.txt", "15Kfitted_fil.txt"),nrow=1)	#the filtered output files
labels<-matrix(data=c('R','2K','15K'),nrow=1)							#suffix to build the output count tables

# .. set the specific input and outputfiles, and call the script that does the cleansing and table building
for (i in 1:3) {
	filename2<-files[i]
	filename_fil<-files_fil[i]
	label<-labels[i]
	source('filter_and_count.R')
}

# .. unload the library to prevent clashes with other functions
detach(package:dplyr)
