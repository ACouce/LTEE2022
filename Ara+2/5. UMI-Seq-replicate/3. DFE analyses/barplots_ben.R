# ..  This script plots the most beneficial alleles in the ancestral background, and their fitness effects in the 2K and 15K backgrounds

# .. load beneficial tails using thben as threshold for approximate neutrality
thben<-0.015
R_beneficial<-Rtable[which(Rtable$fitted1>thben),]
R_beneficial<-R_beneficial[!duplicated(R_beneficial$site),] 			# to impede overlapping regions being counted twice!

# .. preparing objects for a loop to create a table with effects of same allele in Anc, 2K and 15K
RORFs<-as.character(R_beneficial$ORF)
Rlist<-as.character(R_beneficial$alle)
Rtabs<-as.data.frame(matrix(0, length(Rlist), 3))

for (i in 1:length(Rlist)) {
	orf<-R_beneficial[which(R_beneficial$alle==Rlist[i]),]$ORF

	Rtabs[i,1]<-Rtable[which(Rtable$alle==Rlist[i]),]$fitted1
		
	if (length(Ttable[which(Ttable$alle==Rlist[i]),]$fitted1)>0) {		# if exact same allele is available, use it..
		Rtabs[i,2]<-Ttable[which(Ttable$alle==Rlist[i]),]$fitted1	
	} else {Rtabs[i,2]<-mean(Ttable[which(Ttable$ORF==orf),]$fitted1)}	# ..otherwise, take locus average	
		
	if (length(Ftable[which(Ftable$alle==Rlist[i]),]$fitted1)>0) {		# if exact same allele is available, use it..
		Rtabs[i,3]<-Ftable[which(Ftable$alle==Rlist[i]),]$fitted1					
	} else {Rtabs[i,3]<-mean(Ftable[which(Ftable$ORF==orf),]$fitted1)}	# ..otherwise, take locus average
}

# .. add useful columns to the table
Rtabs$name<-RORFs 
Rtabs$site<-R_beneficial$site
Rtabs$abn<-R_beneficial$abn
Rtabs$block<-NA

# .. add column names
colnames(Rtabs) <- c('R', 'M', 'K', 'name', 'site','abn','block')

# .. sort by site (useful for clustering)
Rtabs<-Rtabs[order(Rtabs$site),]

# .. load local context information
locctx<-read.table(file = 'R606_localctx.txt', sep = "\t", header = FALSE)
colnames(locctx) <- c('name', 'str', 'end', 'std', 'ctx')

# .. create a matrix to store block labels (clustered by transcriptional units, as per locctx)
g<-matrix(NA,dim(Rtabs)[1],1)
k<-0

# .. a loop to assign alleles to same block if they belong to the same TU
for (i in 1:dim(Rtabs)[1]) {

	exctname<-paste("\\b",Rtabs[i,]$name,"\\b", sep="")	# trick to ensure perfect match (grep is promiscuous)

	if (!is.na(g[i])) {					# if already assigned a block, keep it and assign it to its neighbours (de facto merging contiguous TUs)

		partners<-which(Rtabs$name%in%locctx[grep(exctname,locctx$ctx),]$name)
		g[partners]<-g[i]
			
	} else {
		if (length(Rtabs[Rtabs$name%in%locctx[grep(exctname,locctx$ctx),]$name,]$name)<=1) { 	# if is a singleton, assign it to a new block
			k<-k+1
			g[i]<-k
		} else {										# otherwise, add their neighbours too to a new block
			partners<-which(Rtabs$name%in%locctx[grep(exctname,locctx$ctx),]$name)
			if (any(!is.na(g[partners]))) {							# to deal with cases in which an allele has been already assigned to a contiguous TU
				g[partners]<-min(g[partners], na.rm = TRUE)
			} else {
				k<-k+1
				g[partners]<-k		
			}
		}	
	}
}

# .. add blocks to table
Rtabs$block<-g

# .. lets sort the table by the fitness of each block 
# .. first, a loop to compute fitness of each block (upper quartile)
blck_fit<-matrix(0, nrow(Rtabs), 1)
for (i in 1:nrow(Rtabs)) {

	blck_fit[i]<-quantile(Rtabs[Rtabs$block%in%(Rtabs[i,]$block),]$R, 0.75)
	maxpos<-which.max(Rtabs[Rtabs$block%in%(Rtabs[i,]$block),]$R)
	if ((Rtabs[Rtabs$block%in%(Rtabs[i,]$block),][maxpos,]$abn)<2) { 	# to discard blocks lead by a low reliability allele	
		Rtabs[Rtabs$block%in%Rtabs[i,],]$R<--100
	}
}

# .. second, blck_fit to do the sorting 
Rtabs$val<-blck_fit
Rtabs<-Rtabs[order(Rtabs$val, decreasing=TRUE),]

# .. eliminate the discarded blocks above
Rtabs<-subset(Rtabs[Rtabs$R>thben,])

# .. rename blocks in order to help with visualization
rename<-matrix(NA, nrow(Rtabs), 1)
r<-0
for (i in 1:nrow(Rtabs)) {
	if (!is.na(rename[i])) {		# if already assigned a block, keep it and assign it to its neighbours
		next }	
	r<-r+1
	rename[which(Rtabs$block%in%(Rtabs[i,]$block))]<-r
}
Rtabs$block<-rename

# .. to mitigate false negatives due to low coverage, let's rescue singletons by using local information (criterion: presence of other initially filtered beneficial mutations within the same block)
singletons<-which(as.vector(table(Rtabs$block))==1)
rescued<-Rtabs[Rtabs$block%in%singletons,]

# .. load unfiltered dataset and apply less stringent filters than initially ("filter_and_count.R") 
RtableU<-read.table(file = "Rfitted_log.txt", sep="\t", header=TRUE, as.is=TRUE, comment.char = "") # to change accordintly to obtain the lists of beneficial mutations for 2K and 15K needed in folder "4. Predictability"
RtableU<-RtableU[!is.na(RtableU$sterr1),]
RtableU[(RtableU$t2<=5 | RtableU$sterr1>0.02) & RtableU$fitted1>0.015,]$sterr1 <- NA
RtableU<-RtableU[!is.na(RtableU$sterr1),]

rlist<-c(NA,NA, NA)
for (r in 1:nrow(rescued)) {

	nei<-unlist(strsplit(as.character(locctx[locctx$name==rescued[r,]$name,]$ctx), split=","))
	if (sum(RtableU[RtableU$ORF%in%nei,]$fitted1>0.015)-1 >= 1) {				# rescue if there is at least another allele supporting the beneficial effect
		neis<-RtableU[RtableU$ORF%in%nei & RtableU$fitted1>0.015,]
		for (u in unique(neis$alle)) {
			saved<-c(u,rescued[r,]$block, rescued[r,]$val)
			rlist<-rbind(rlist,saved)
		}
	}
}
rlist<-as.data.frame(rlist[-1,])

# .. recycle RtabsU to store effects of rescued allele in Anc, 2K and 15K
RtabsU<-as.data.frame(matrix(0, nrow(rlist), 8))
for (i in 1:nrow(rlist)) {

	orf<-RtableU[which(RtableU$alle==rlist[i,1]),]$ORF
	
	RtabsU[i,1]<-RtableU[which(RtableU$alle==rlist[i,1]),]$fitted1

	if (length(Ttable[which(Ttable$alle==rlist[i,1]),]$fitted1)>0) {			# if exact same allele is available, use it..
		RtabsU[i,2]<-Ttable[which(Ttable$alle==rlist[i,1]),]$fitted1
	
	} else {RtabsU[i,2]<-mean(Ttable[which(Ttable$ORF==orf),]$fitted1)}		# ..otherwise, take locus average
	
	if (length(Ftable[which(Ftable$alle==rlist[i,1]),]$fitted1)>0) {			# if exact same allele is available, use it..
		RtabsU[i,3]<-Ftable[which(Ftable$alle==rlist[i,1]),]$fitted1
				
	} else {RtabsU[i,3]<-mean(Ftable[which(Ftable$ORF==orf),]$fitted1)}		# ..otherwise, take locus average
		
	RtabsU[i,4]<-orf
	RtabsU[i,5]<-RtableU[which(RtableU$alle==rlist[i,1]),]$site
	RtabsU[i,6]<-RtableU[which(RtableU$alle==rlist[i,1]),]$abn
	RtabsU[i,7]<-as.numeric(as.character(rlist[i,2])) 				# to avoid issues with factors
	RtabsU[i,8]<-as.numeric(as.character(rlist[i,3]))

}
colnames(RtabsU) <- c('R', 'M', 'K', 'name', 'site','abn','block','val')

# .. prepare final table to plot by merging the non singletons (i.e., "multiple") with the rescued singletons
multiple<-which(as.vector(table(Rtabs$block))>1)
toplot<-rbind(Rtabs[Rtabs$block%in%multiple,], RtabsU)

# .. let's only show the average fitness per locus (instead of individual values for each segment)
benorf<-unique(toplot$name)
tsingle<-as.data.frame(matrix(0, length(benorf), 8))
n<-0
for (i in benorf) {
		n<-n+1
		tsingle[n,1:3]<-apply(toplot[toplot$name%in%i,1:3], 2, mean,na.rm=TRUE)	# average fitness per locus 
		tsingle[n,4:8]<-toplot[toplot$name%in%i,4:8][1,]				# just copy details
	}
toplot<-tsingle[,1:8]

# .. add column names
colnames(toplot) <- c('R', 'M', 'K', 'name', 'site','abn','block','bvalue')

# .. to update block labels (we missed some along the way)
oldblocks<-unique(toplot$block)
newblocks<-1:length(oldblocks)
for (b in 1:length(oldblocks)) {
	toplot[toplot$block==oldblocks[b],]$block<-newblocks[b]
}

# .. store data table for later analyses ('metagenomics.R')
write.table(toplot,'Rtoanal.dat', sep = "\t") 	# to obtain the lists of beneficial mutations needed in folder "4. Predictability", change Rtable for Ttable and Ftable, load RtableU accordingly, and relabel here

# .. initialize graphical output
png("bens1.png", width=1100*0.011, height=(500*0.865)*0.011, units = 'in', res = 300)
par(family="sans", cex=1)
par(mar = c(8.1, 4.1, 4.1, 2.1)) 

# .. show the most beneficial 55 loci in the ancestor (22 blocks, cut just for visualization purposes)
toplot<-toplot[toplot$block<toplot[55,]$block,]

# .. let's make a barplot and store coordinates to add white and gray shaded areas grouping loci from same block
coord<-barplot(t(toplot[,1:3]), beside = TRUE, space=c(0,1), axes = FALSE, axisnames=FALSE, ylim=c( -0.08, 0.11), col='white', las=2, cex.names=1, cex.axis=1.2, horiz = FALSE, add = FALSE, border = 'NA', lwd=1.2)
coord<-rbind(coord,t(toplot$block))

# .. adding white and gray shaded areas
for (r in seq(1, max(toplot$block), by=2)) {
	rect(min(coord[1:3,coord[4,]==r+1])-1,-3,max(coord[1:3,coord[4,]==r+1])+1,3,col="white",lty=1, lwd=0.5)
	rect(min(coord[1:3,coord[4,]==r])-1,-3,max(coord[1:3,coord[4,]==r])+1,3,col="grey85",lty=1, lwd=0.5)
	}

# .. the actual barplot
barplot(t(toplot[,1:3]), beside = TRUE, space=c(0,1), yaxt='n', names.arg=toplot[,4], font=3, ylim=c(-0.08, 0.11), col=c('black', col2K, col15K), las=2, cex.names=1.1, cex.axis=1.2, horiz = FALSE, add = TRUE, lwd=1)

# ..  create matrix to store average standard errors
ster<-matrix(NA,nrow(toplot),3)
for (s in 1:nrow(toplot)) {
	name<-toplot[s,]$name
	ster[s,1]<-mean(Rtable[Rtable$ORF==name,]$sterr1, na.rm=TRUE)
	ster[s,2]<-mean(Ttable[Ttable$ORF==name,]$sterr1, na.rm=TRUE)
	ster[s,3]<-mean(Ftable[Ftable$ORF==name,]$sterr1, na.rm=TRUE)
}

# .. use standard errors to add errors bars indicate ~95% confidence intervals
segments(coord[1,], toplot$R-2*ster[,1], coord[1,], toplot$R+2*ster[,1], lwd=1.2)
segments(coord[2,], toplot$M-2*ster[,2], coord[2,], toplot$M+2*ster[,2], col=col2K, lwd=1.2)
segments(coord[3,], toplot$K-2*ster[,3], coord[3,], toplot$K+2*ster[,3], col=col15K, lwd=1.2)

# .. custom y-axis
axis(2, at = seq(-0.075, 0.1, 0.025), labels=NA, lwd=1.5, cex.axis = 1, tck=-0.05)
axis(2, at = seq(-0.05, 0.1, 0.05),  labels=seq(-0.05, 0.1, 0.05), lwd=1.5, cex.axis=1.35, tck=-0.05)
abline(h=0, lty=2)

# .. close graphical output
dev.off()

# .. initialize graphical output
png("bens2.png", width=300*0.8*0.011, height=460*0.85*0.011, units = 'in', res = 300)
par(family="sans", cex=1)

# .. plot summary barplots
boxplot(toplot[toplot$R>0.015,1:3], col=c('grey24', col2K, col15K), ylim=c(-0.08, 0.11), cex.names=1, cex=0.8, names=c('Anc', '2K', '15K'), ylab='selection coeff. (s)', cex.lab=1.5,  cex.axis=1.2, lwd=1, yaxt='n', las=2)
abline(h=0, lty=2)

# .. custom y-axis
axis(2, at = seq(-0.075, 0.1, 0.025), labels=NA, lwd=1.5, cex.axis = 1, tck=-0.05)
axis(2, at = seq(-0.05, 0.1, 0.05),  labels=seq(-0.05, 0.1, 0.05), lwd=1.5, cex.axis=1.35, tck=-0.05)

# .. homogenize border lines
box(lwd=1.5)

# .. close graphical output
dev.off()
