# .. This script filters out low-quality fits and presumed hitchhikers, and then constructs the count tables required for plotting the DFE.

# .. define a useful function
'%!in%' <- function(x,y)!('%in%'(x,y))

# .. load the specific dataset to work with
ttable<-read.table(file = filename2, sep="\t", header=TRUE, as.is=TRUE, comment.char = "")

# .. remove NAs, alleles with very low initial abundance, and beneficial alleles with large regression standard errors
ttable[is.na(ttable$sterr1),c(4:24)]<--100
ttable[ttable$t1<10,c(13:24)]<--101
ttable<-subset(ttable[ttable$fitted1>-1,])
ttable[ttable$sterr1>0.01 & ttable$fitted1>0,c(13:24)]<--102
filter<-subset(ttable[ttable$fitted1>-1,])

# .. preparing objects for a loop to eliminate presumed hitchhikers by comparing fitness effects across segments of the same locus
# .. first, singletons are removed, as they cannot be filtered by comparison, and only high-quality variants are retained (those with at least 4 independent observations and initial abundance counts exceeding 20)
tab<-table(filter$ORF)
reps<-names(tab[tab>1])
singletons<-names(tab[tab==1])
filter[filter$ORF%in%singletons & filter$abn<=3 & filter$t1<=20 & filter$fitted1>0,13:24]<--107

# .. create reference list of genes to monitor progress (clockwise from origin)
linkage<-c('lacZ-1', 'galE-1', 'trpA-1', 'hipA-1', 'cheA-1', 'gyrA-1', 'purL-1', 'argR-1', 'xylA-1', 'uvrA-1') 

for (r in  reps) {	# loop accross the filtered gene names (skipping loci for which functional similarity accross internal fragments is not guaranteed, i.e., intergenic regions, pseudogenes, etc.)

	if (length(grep('IR',r))==1) {next}	
	if (length(grep('ESCRE',r))==1) {next}
	if (length(grep('ECB',r))==1) {next}
	if (length(grep('#',r))==1) {next}	

	alls<-filter[filter$ORF%in%r,]				# load all alleles from same loci
	if (alls$alle[1]%in%linkage) {print(r)} 			# report loop progress

	if (length(grep(r,alls$alle)) != length(alls$alle)) {	# make sure to consider only segments of the same locus (discarding gene overlaps)
		alls<-alls[grep(r,alls$alle),]
		if (nrow(alls)<2) {next} 			# discard gene overlaps
	}
	
	dev<-as.numeric(which.max(colSums(as.matrix(dist(alls$fitted1, upper=TRUE)))))	# identify most allele most deviated within locus
	ave<-mean(alls$fitted1)								# estimate locus average
	
	if (abs(alls[dev,]$fitted1-ave)/abs(mean(alls$fitted1))>=0.25) { 			# if deviation represents > 25% of locus average

		alls_wo<-alls[!(grepl('-1',alls$alle) & alls$strand[1]=='C') & !(grepl('-5',alls$alle) & alls$strand[1]=='F'),]	# first, check whether deviation is due to N- or C-terminal segments
		dev<-as.numeric(which.max(colSums(as.matrix(dist(alls_wo$fitted1, upper=TRUE)))))					# recalculate most deviated allele without terminal segments 
		ave<-mean(alls_wo$fitted1)											# recalculate locus average without terminal segments 
		if (abs(alls_wo[dev,]$fitted1-ave)/abs(mean(alls$fitted1))<=0.25) { 						# if deviation now falls < 25% of locus average
				
			ter<-alls[which((grepl('-1',alls$alle) & alls$strand[1]=='C') | (grepl('-5',alls$alle) & alls$strand[1]=='F')),]$alle
			if (length(ter)>0) {
				filter[filter$alle%in%ter & filter$abn<=3 & filter$t1<=20,13:24]<--107 				# retain the terminal segment only if reliable enough			
			}
			next													# no further outliers, move to next gene name	
		}
				
		if (nrow(alls)==2) {	# in cases where only two segments are present, keep the most reliable based on number of independent observations ($abn) and initial abundance counts ($t1)
			
			aux<-alls[alls$abn>3 & alls$t1>20,c(2,5)]
			chosen<-aux[which.max(aux$t1),]$alle			# if both are good, keep the one with higher initial abundance counts
			if (length(chosen)==0) {
				filter[filter$ORF%in%r,13:24]<--103
			} else {		
				filter[filter$alle%in%alls[alls$alle%!in%chosen,]$alle,13:24]<--103
			}
			next							# no further outliers, move to next gene name	
		}
		
		outs<-alls[!between(alls$fitted1, ave-0.01, ave+0.01),]$alle	# indentify alleles that deviate by > 0.01 in fitness from the locus average
		
		ter<-outs[which((grepl('-1',outs) & alls$strand[1]=='C') | (grepl('-5',outs) & alls$strand[1]=='F'))]	# check whether deviation is due to N- or C-terminal segments
		if (length(ter)>0) {
			filter[filter$alle%in%ter & filter$abn<=3 & filter$t1<=20,13:24]<--107 				# retain the terminal segment only if reliable enough
			outs<-outs[outs!=ter]
		}

		filter[filter$alle%in%outs,13:24]<--105				# remove the deviated allele
	}
}

# .. update data table
ttable<-filter

# .. set lower limit for deleterious effects
ttable[ttable$fitted1<=-0.5,]$fitted1<-NA
	
# .. store clean data table for further reference
write.table(ttable, file = filename_fil, col.names = TRUE, sep="\t")

# .. create objects for constructing a count table of fitness classes (required for plotting the DFE)
r_fitted_count<-matrix(0,nrow = 60, ncol = 4)
orgn=-0.42
step=0.01
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
